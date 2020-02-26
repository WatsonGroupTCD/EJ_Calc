!This program requires a file having the spin densities for the system of
!interest in different electronic configurations and the energy associated with
!the system in that electronic configuration. This file can be generated using
!the program spin_dens_extrctr. Besides this, the program requires another file
!having details about the number of magnetic centres and the Hamiltonian. Using
!these two files, this program lets one calculate all possible J values for the
!system. If you have data for n configurations and you need j coupling constants
!then this program calculates all nCj coupling constants (provided n>j). In each
!calculation you need to sacrifice one equation. Thus you basically end up with
!n sets of (n-1)Cj coupling constants, all of which are calculated. The average
!per set and the total average of the coupling constant is also calculated.

Module modul_lib
 Use init
 Use enrg_frm_jval
 Use Iso_Fortran_Env
 Implicit None
        
 Contains

 !This subroutine will determine the S1S2 terms and how many of each term comes
 !under 1 J value
 Subroutine jvals()
        Integer(kind = int_kind) :: i1, num_j_chk
        Allocate (jval(num_spin_dens_set, Size(jpos) - 1), tot_avg(Size(jpos) - 1)) 
        j = 2
   !!!!! Initialising jval matrix to zero !!!!!
        Do i1 = 1, num_spin_dens_set
                Do i = 1, Size(jpos) - 1
                        jval(i1, i) = 0
                End Do
        End Do
   !!!!! Initialization Ends !!!!!
 
   !!!!! Determining the coefficient of J values !!!!!
        Do i1 = 1, num_spin_dens_set
                j = 2
                Do i = 1, tot_interac
                        If (hamil1(i) .ne. 0 .AND. hamil2(i) .ne. 0) then
                        If (i == jpos(j)) Then
                                j = j + 1
                        End If 
                        jval(i1, j - 1) = jval(i1, j - 1) + hamil(i1,i)
                        End If
                End Do
        End Do
        num_j_chk = 0
        Do i = 1, num_spin_dens_set
                If (jval(i, Size(jpos) - 1) == 0) then
                        num_j_chk = num_j_chk + 1
                End If
        End Do
        If (num_j_chk == num_spin_dens_set) then
                Write(*,*)'You have specified more J values but have defined less'
                Write(*,*)'Aborting'
                Call abort()
        End If

        Write(12, *) 'Coefficients of J values' 
        178 Format(F10.4)
        155 Format(i2)
        Do i = 1, Size(jpos) - 1
                Write(12, '(a)', advance = 'no') '    J'                
                Write(12, 155, advance = 'no')i
                Write(12, '(a)', advance = 'no') '   '
        End Do
        Write(12, *)
        Do i1 = 1, num_spin_dens_set
                Do i = 1, Size(jpos) - 1
                        Write(12, 178, advance = 'no') jval(i1, i) 
                End Do
                Write(12, *) 
        End Do
 
 End Subroutine jvals
 !The following subroutine is the main routine. It forms all possible sets of
 !equations and solves them 
 Subroutine solv()
        Integer(kind = int_kind) :: k, l, m, counter, x, z
        Integer(kind = int_kind), Dimension(:), Allocatable :: ctr
        !The matrix ctr will allow the seperation and hence the evaluation of
        !all possible combinations.
        Real(kind = real_kind), Dimension(:), Allocatable :: loc_energy, enrg_tb_slvd, enrg_tb_slvd1
        Real(kind = real_kind), Dimension(:,:), Allocatable :: loc_jval, jval_tb_solvd, jval_tb_solvd1        
        !loc_energy and loc_jval are used to locally store the coefficients of j
        !values and energy since they will be changed in each set of equations
        !considered.
        !jval_tb_solvd1 and enrg_tb_slvd1 are essentially copies of the arrays
        !jval and energy after deducting the reference equation from all of them
        !but the refernce equation that will be used has been removed from this
        !set.

        !jval_tb_solvd and enrg_tb_slvd store the
        !set of equations that will be solved in one instance

        Logical :: chk
        Allocate(loc_jval(num_spin_dens_set, Size(jpos) - 1), loc_energy(num_spin_dens_set))
        Allocate(ctr(Size(jpos)))
        Allocate(jval_tb_solvd(Size(jpos) - 1,Size(jpos) - 1), enrg_tb_slvd(Size(jpos) - 1))
        Allocate(jval_tb_solvd1(num_spin_dens_set - 1,Size(jpos) - 1), enrg_tb_slvd1(num_spin_dens_set - 1))
        x = num_spin_dens_set - 1
        poss_comb = comb(x,(Size(jpos) - 1))
        Write(12,*) 'Possible combination of equations per set of equations '
        Write(12,*) 'considered (i.e. total number of eqations -1):' 
        Write(12,*) x, 'C', (Size(jpos) - 1), '=', poss_comb
        Allocate(jval_trkr(num_spin_dens_set, poss_comb, Size(jpos) - 1))
        Allocate (singul(num_spin_dens_set, poss_comb))
        789 Format(i5) 
        singul = 0
        Do i = 1, num_spin_dens_set
                Read (14,*) (singul(i,j),j=1, poss_comb)
        End Do
        Write (12,*) 'Solving the equations'
        Do i = 1, num_spin_dens_set ! This loop changes the reference equation
                Do z = 1, Size(jpos) - 1
                        ctr(z) = z
                End Do
                ctr(z) = num_spin_dens_set
                chk = .False.
                loc_energy = energy
                loc_jval = jval
                loc_energy = (loc_energy - energy(i)) * 2.19474625 * 10**5 
                !The factor *627.5095*4.18*83.53548 was previously being used
                !for the conversion from hartree to cm-1


                !The coefficients of the reference equation are being deducted
                !from all the equations present in the following loop section 
                Do k = 1, num_spin_dens_set
                        Do l = 1, Size(jpos) - 1
                                loc_jval(k,l) = loc_jval(k,l) - jval(i,l)
                        End Do
                End Do

                !The following loop section removes the reference equations
                !(which has all the coefficients and energy equated to zero)
                !from the set of all possible equations. The new set formed will
                !be used for further calculations.
                counter = 1
                Do k = 1, num_spin_dens_set
                        Do l = 1, Size(jpos) - 1                
                                if (loc_energy(k) .ne. 0) then
                                        chk = .True.
                                End If
                        End Do
                        If (chk == .True.) Then
                                Do l = 1, Size(jpos) - 1
                                        jval_tb_solvd1(counter,l) = loc_jval(k,l)
                                End Do
                                enrg_tb_slvd1(counter) = loc_energy(k)
                                counter = counter + 1
                        End If
                        chk = .False.
                End Do
   !!!!! Printing Area !!!!!
                112 Format(2x, F10.3)
                !Printing the new set of equation wrt the current reference
                !state
                !Uncomment the following only for debugging. This will give you
                !the coefficient of J values with respect to the reference
                !equation.
                !Write (12,*)
                !Write (12,*) 'Weight Coefficients for the reference state:'
                !Write (12, 112, advance = 'no') (jval(i,l), l = 1, Size(jpos) - 1)
                !Write (12,*)
                !Write (12,*) 'Energy of reference state:', i
                !Write(12, 112, advance= 'no') energy(i)
                !Write (12,*)
                !Write (12,*) 'Coefficients and Energy wrt the reference state'
                !Do k = 1, num_spin_dens_set
                !        Do l = 1, Size(jpos) - 1
                !                Write(12, 112, advance= 'no') loc_jval(k,l)
                !        End Do
                !        Write(12, 112, advance= 'no') loc_energy(k)
                !        Write(12,*)
                !End Do
                !!Printing the same set of equations with the current reference
                !!state removed since that will be zero
                !Write (12, *) 'The set of equations for which all possible'
                !Write (12, *) 'combinations will be considered: '
                !Do k = 1, num_spin_dens_set - 1 
                !        Do l = 1, Size(jpos) - 1
                !                Write (12, 112, advance='no') jval_tb_solvd1(k,l)
                !        End Do
                !        Write(12, 112, advance= 'no') enrg_tb_slvd1(k)
                !        Write(12,*)
                !End Do
                !Write (12,*)
    !!!!! Printing Ends !!!!!
                !The following section will go through all possible combinations
                !and solve each of them to obtain the J-values.
                counter = 1 
                !counter keeps track of the number of rows
                !jval_tb_solvd and enrg_tb_slvd matrices will have. These
                !matrices contain the set of equations that are solved
                !simultaneously. If you need 2 j values then these matrices will
                !have two rows each since that is the number of equations
                !required to determine 2 unknowns.jval_tb_solvd will be a 2X2
                !matrix since it will have coefficients for both j
                !terms. enrg_tb_slvd matrix just gives the energy for a
                !particular configuration.
                x = num_spin_dens_set - 1 !x is the total number of non-zero
                                          !equations which is obviously 1 less
                                          !than the total equations
                Do k = 1, poss_comb
                !i.e. the total possible combinations of equations
                        Do m = 1, Size(jpos) - 1 
                                Do l = 1, Size(jpos) - 1
                                        jval_tb_solvd(counter, l) = jval_tb_solvd1(ctr(m), l)
                                        enrg_tb_slvd(counter) = enrg_tb_slvd1(ctr(m))
                                End Do        
                                        counter = counter + 1
                                If (counter == Size(jpos)) then
                                        counter = 1
                                End If
                        End Do
   !!!!! Printing Area !!!!!                               
                        !Debugger. Prints the set of equations that are being
                        !solved.

                        !If(singul(i,k) == k) then
                        !Write (12,*)k, '- j val set which will be solved'
                        !Do z = 1, Size(jpos) - 1 
                        !        Do l = 1, Size(jpos) - 1
                        !                Write (12, 112, advance='no') jval_tb_solvd(z,l)
                        !        End Do
                        !        Write(12, 112, advance= 'no') enrg_tb_slvd(z)
                        !        Write(12,*)
                        !End Do
                        !End If
   !!!!! Printing Ends !!!!!
                      If(singul(i,k) == k) then
                        !Write (12,*) 'J values(cm-1) are'
                        Call simul_solv(jval_tb_solvd, enrg_tb_slvd)
                        jval_trkr(i, k, :) = enrg_tb_slvd
                      End If
                        ctr(Size(ctr)-1) = ctr(Size(ctr)-1) + 1
                        Call checkpos(ctr)
                        If (counter == Size(jpos)) then
                                counter = 1
                        End If
                End Do
        End Do   
 EndSubroutine solv

 !The checkpos subroutine takes care of the fact that ctr cycles through only
 !the possible combinations and does not go out of bounds.
 Subroutine checkpos(ctr1)
        Integer(kind = int_kind), Dimension(:), Allocatable :: ctr1
        Integer(kind = int_kind) :: k, l, m, max_at_pos
        If (ctr1(Size(ctr1)-1) == ctr1(Size(ctr1))) then
                Do k = Size(ctr1)-2, 1, -1
                        max_at_pos = ctr1(Size(ctr1)) - (Size(ctr1) - k)
                        If (ctr1(k) <= max_at_pos - 1) then
                                ctr1(k) = ctr1(k) + 1
                                Do l = k+1, Size(ctr1)-1
                                        ctr1(l) = ctr1(l-1) + 1
                                End Do
                        Exit
                        End If
                End Do
        End If
 End Subroutine checkpos
 
 !The following subroutine solved the set of equations generated by the solv
 !subroutine and solves them using a LAPACK routine
 Subroutine simul_solv(jval_t_solv, enrg_t_slv)
        Integer(kind = int_kind) :: n, lda, info, nrhs, ldb
        Integer(kind = int_kind), Dimension(:), Allocatable :: ipiv
        Real(kind = real_kind), Dimension(:), Allocatable :: enrg_t_slv 
        Real(kind = real_kind), Dimension(:,:), Allocatable :: jval_t_solv
        n = Size(enrg_t_slv)
        lda = Size(enrg_t_slv)
        nrhs = 1
        ldb = n
        Allocate(ipiv(n))
        !look at the following website
        !http://physics.oregonstate.edu/~landaur/nacphy/lapack/fortran.html
        Call dgesv (n, nrhs, jval_t_solv, lda, ipiv, enrg_t_slv, ldb, info)
        !The following will print the J value set obtained by solving the given
        !set of equations.
        !Do n = 1, Size(enrg_t_slv)
        !        134 Format(i1)
        !        112 Format(2x, F18.10)
        !        If (Abs(enrg_t_slv(n)) > 7) then
        !                Write (12, *) ' This value looks very big!!!'
        !        End If
        !        Write (12,'(a)', advance = 'no') 'J'
        !        Write (12, 134, advance = 'no') n
        !        Write (12,'(a)', advance = 'no') '='
        !        Write (12, *) enrg_t_slv(n)
        !End Do
        !Write (12,*) 'Info', info
        !Write (12,*)
 End Subroutine simul_solv
 !The following subroutine calculated the average of all the valid set of
 !solutions
 Subroutine avrg_calcs()
        Integer(kind = int_kind) :: k, l, m, x, poss_comb, counter
        Real(kind = real_kind), Dimension(:), Allocatable :: avg
        Allocate(avg(Size(jpos) - 1))
        counter = 0
        main_cntr = 0
        avg = 0
        tot_avg = 0
        112 Format(1x, F10.5)
        134 Format(i2)
        Write (12,*) ' Average J values'
        Do m = 1, Size(jpos) - 1
                Write (12,'(a)', advance = 'no') '     '
                Write (12,'(a)', advance = 'no') 'J'
                Write (12, 134, advance = 'no') m
                Write (12,'(a)', advance = 'no') '    '
        End Do
        Write (12,*)
        Write (12,*) ' Average across each reference equation set:'
        x = num_spin_dens_set - 1
        poss_comb = comb(x,(Size(jpos) - 1)) 
        Do k = 1, num_spin_dens_set
                Do l = 1, poss_comb
                        If (jval_trkr(k, l,1) .ne. 0) then
                        counter = counter + 1
                        End If
                        Do m = 1, Size(jpos) - 1
                                If(jval_trkr(k, l, m) .ne. 0) Then
                                        avg(m) = avg(m) + jval_trkr(k, l, m) 
                                        !counter = counter + 1
                                End If
                        End Do
                End Do
                Write(12,*) 'Total non singular solutions = ', counter
                avg = avg/(counter) !poss_comb
                main_cntr = main_cntr + counter
                counter = 0
                Do m = 1, Size(jpos) - 1
                        Write(12, 112, advance = 'no') avg(m)
                End Do
                Write (12,*)
                tot_avg = tot_avg + avg
                avg = 0
        End Do
        tot_avg = tot_avg/num_spin_dens_set
        Write(12,*) ' Global Average'
        Do m = 1, Size(jpos) - 1
                Write(12,112, advance = 'no') tot_avg(m)
        End Do
        Write (12,*)
 End Subroutine avrg_calcs

 !!!Implementation of standard deviation!!!
         !This will calculate the standard deviation for each J value
         Subroutine std_dev(avg)
                Integer(kind = int_kind) :: k, l, m
                Real(kind = real_kind):: ss
                Real(kind = real_kind), Dimension(:), Allocatable :: avg, var, stddev, std_per
                Allocate (var(Size(avg)), stddev(Size(avg)), std_per(Size(avg)))
                var = 0
                stddev = 0
                Do i = 1, num_spin_dens_set
                        Do j = 1, poss_comb
                                Do k = 1, Size(avg)
                                        If (jval_trkr(i, j, k) .ne. 0) then
                                                var(k) = var(k) + (avg(k) - jval_trkr(i, j, k))**2
                                        End If
                                End Do
                        End Do
                End Do
                var = var/main_cntr
                stddev(:) = sqrt(var(:))
                std_per(:) = stddev(:) / Abs(avg(:)) * 100
                !stdv = .TRUE
                Write(12,*)'total non-singular equations', main_cntr
                Write(12,*)'standard deviation (percentage standard deviation):'
                165 Format(1x, F6.2)
                166 Format(i3)
                Do k = 1, Size(avg)
                        Write(12, 165, advance = 'no') stddev(k)
                        Write(12, '(a)', advance = 'no') '('
                        Write(12, 166, advance = 'no') Int(std_per(k))
                        Write(12, '(a)', advance = 'no') '%)'
                End Do
                Write(12,*)
                Call std_dev1(stddev)
         End Subroutine std_dev
         !This will check if all the valid J value sets are within 3 standard
         !deviations or not. If they are not, they will be eliminated
         Subroutine std_dev1(stdev)
                Integer(kind = int_kind) :: k, l, m, tmp_cntr
                Real(kind = real_kind), Dimension(:), Allocatable :: stdev
                tmp_cntr = main_cntr
                Do i = 1, num_spin_dens_set
                        Do j = 1, poss_comb
                                Do k = 1, Size(stdev)
                                        If (Abs(jval_trkr(i, j, k)) > 0) then
                                                If(jval_trkr(i, j, k) > (tot_avg(k) + 3.0*stdev(k)) .OR. jval_trkr(i, j, k) < (tot_avg(k) - 3.0*stdev(k))) then
                                                        jval_trkr(i, j, :) = 0.0
                                                        tmp_cntr = tmp_cntr - 1 
                                                End If
                                        End If
                                End Do
                        End Do
                End Do
                If (tmp_cntr == main_cntr) then
                        stdv = .TRUE.
                Else
                        main_cntr = tmp_cntr
                        Call avrg_calcs1()
                End If
         End Subroutine std_dev1
         !This will calculate the new average for J values since some have been
         !eliminated by the above subroutine
         Subroutine avrg_calcs1()
                Integer(kind = int_kind) :: k, l, m, tmp_cntr
                tot_avg = 0
                Do i = 1, num_spin_dens_set
                        Do j = 1, poss_comb
                                Do k = 1, Size(tot_avg)
                                        tot_avg(k) = tot_avg(k) + jval_trkr(i, j, k)
                                End Do
                        End Do
                End Do
                tot_avg = tot_avg/main_cntr
         End Subroutine avrg_calcs1


 Subroutine final_print()
        Integer(kind = int_kind) :: k
        165 Format(1x, F10.5)
        Write(12,*)
        Write(12,*) 'Total non singular equations', main_cntr
        Write(12,*) 'new average'
        Do k = 1,  Size(tot_avg)
                Write(12, 165, advance = 'no') tot_avg(k)
        End Do
        !Debugger. Prints all the J value sets that have been considered valid
        !after applying the standard deviation constraint.

        !Write(12,*) 
        !Write(12,*) 'J value sets'
        !Do i = 1, num_spin_dens_set
        !       Do j = 1, poss_comb 
        !               Do k = 1, Size(tot_avg)
        !                       If (Abs(jval_trkr(i, j, k)) > 0) then
        !                               If( jval_trkr(i, j, k) > 10) then
        !                                    Write(12,*) 'This values is still big', i, j, k 
        !                               End If
        !                               Write(12, 165, advance = 'no') jval_trkr(i, j, k)
        !                       End If
        !               End Do
        !               If (Abs(jval_trkr(i, j, 1)) > 0) then
        !               Write(12,*)
        !               End If
        !       End Do
        !End Do
 End Subroutine final_print

 !This subroutine will calculate the energies using the computed J-values which
 !can be compared with the DFT energies for those equations 
 Subroutine backtrack(coeff, ener)
        Integer(kind = int_kind) :: k, l, m
        Real(kind = real_kind), Dimension(:), Allocatable :: ener, ener1
        Real(kind = real_kind), Dimension(:,:), Allocatable :: coeff
        k = Size(coeff)/Size(ener)
        Allocate (ener1(k))
        Write (12,*)k, '- j val set which will be solved'
        112 Format(1x, F10.5)
        !Do l = 1, k 
        !        Do m = 1, Size(ener)
        !                Write (12, 112, advance='no') coeff(l,m)
        !        End Do
        !        !Write(12, 112, advance= 'no') ener(l)
        !        Write(12,*) 
        !End Do
        ener1 = 0.0
        Do l = 1, k 
                Do m = 1, Size(ener)
                        ener1(l) = ener1(l) + coeff(l,m)*ener(m)
                End Do
        End Do
        Write(12,*) "backtrack"
        !Write(12,*) ener1
        Do l = 1, k
                Write(12,*) ener1(l) - ener1(1)
        End Do
 End Subroutine backtrack
 !This subroutine will calculate the energies using the computed J-values and
 !compares them with the DFT energies (The final part of the output)
 Subroutine backtrack1(coeff, ener)
        Integer(kind = int_kind) :: k, l, m
        Real(kind = real_kind) :: temp_ener, temp_ener1
        Real(kind = real_kind), Dimension(:), Allocatable :: ener, ener1
        Real(kind = real_kind), Dimension(:,:), Allocatable :: coeff
        k = Size(coeff)/Size(ener)
        Allocate (ener1(k))
        113 Format(3x, F10.5)
        112 Format(14x, F10.5)
        134 Format(F4.1)
        ener1 = 0.0
        Do l = 1, k
                Do m = 1, Size(ener)
                        ener1(l) = ener1(l) + coeff(l,m)*ener(m)
                End Do
        End Do
        Write(12,*)
        Write(12,*) 'Comparison of energy (cm-1)  of electronic states calculated by DFT'
        Write(12,*) 'and energy obtained using the avergae J-values calculated:'
        Write(12,*) '(Note: The following are the energies assuming the first'
        Write(12,*) 'energy value in the spin density file as the reference value)'
        Write(12,*) 'DFT Energy(cm-1)    Calculated energy(cm-1)     Absolute Difference(cm-1)'
        Do l = 1, k
                temp_ener = (energy(l) - energy(1)) * 2.19474625 * 10**5
                temp_ener1 = ener1(l) - ener1(1)
                Write(12, 113, advance = 'no') temp_ener 
                Write(12, 112, advance = 'no') temp_ener1
                Write(12, 112, advance = 'no') Abs(temp_ener - temp_ener1)
                Write(12, '(a)', advance = 'no') '('
                Write(12, 134, advance = 'no') 100.0*Abs((temp_ener - temp_ener1)/temp_ener)
                Write(12, '(a)', advance = 'no') '%)'

                Write(12,*)
        End Do
        Write(12,*)
 End Subroutine backtrack1

!Functions!
         Integer Function comb(n,r)
                Integer(kind = int_kind) :: n, r, k, cntr
                comb = 1
                If (n-r > r) then
                        cntr = r
                        Do k = n-r+1, n
                                comb = comb * k
                                If(comb > 10000 .AND. cntr > 1) then
                                        comb = comb/cntr
                                        cntr = cntr - 1
                                End If
                        End Do
                        comb = comb / fact(cntr)
                Else
                        cntr = n-r
                        Do k = r+1, n
                                comb = comb * k
                                If(comb > 10000 .AND. cntr > 1) then
                                        comb = comb/cntr
                                        cntr = cntr - 1
                                End If
                        End Do
                        comb = comb / fact(cntr)
                End If
                Return
         End function comb
         !Function to calculate the factorial
         Integer Function fact(facto)
                Integer(kind = int_kind) :: facto, k, l
                l = facto
                fact = 1
                Do k = l, 1, -1
                        fact  = fact * k
                End Do
                return
         End function fact

End Module modul_lib 
