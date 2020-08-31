Module init

 Use Iso_Fortran_Env
 Implicit None
 Integer, parameter :: int_kind = INT32
 Integer, parameter :: real_kind = REAL64 !Double precision
 Integer, parameter :: comp_kind = REAL64        
 Integer(kind = int_kind) :: i, j, num_mag_cent, no_of_j_val, tot_interac
 !i and j are for running loops
 !num_mag_cent stores the number of metal centres in each configuration
 !no_of_j_val stores the total number of Jvalues defined for the system
 !tot_interac stores the number of possible interactions between spins i.e. the
 !total number of S1S2 terms.
 Integer(kind = int_kind) :: ios, num_spin_dens_set, try, try1, poss_comb
 !num_spin_dens_set stores the number of sets of spin density available
 !which is the total numberof lines in the spin density file
 !poss_comb determined the possible combination of equations i.e. the
 !(n-1)Cj term (read line 10 from top)
 !try and try1 help in keeping the program running just in case there is an
 !equation set is encountered that cannot be solved.
 Integer(kind = int_kind), Dimension(:), Allocatable :: hamil1, hamil2, jpos
 !jpos keeps track of how many interactions are being included under a single J
 !value
 !hamil1 and hamil2 together keep track of which centres are interacting
 Integer(kind = int_kind) :: main_cntr
 !main_cntr will keep track of total non singular equations
 Real(kind = real_kind), Dimension(:), Allocatable :: energy, tot_avg
 !energy stores the energy for each electronic configuration
 !tot_avg stores the average of all the J-value sets which are valid and which
 !meet the standard deviation criteria.
 Real(kind = real_kind), Dimension(:), Allocatable :: stddev, std_per
 !stddev and std_per store the standard deviation and percentage standard
 !deviation for the different J-values
 Real(kind = real_kind), Dimension(:,:), Allocatable :: spin, hamil, jval
 !spin stores the spin density value associated with each metal centre as
 !given in the spin density file.
 !hamil stores the 2S1S2 type weight terms for each interaction. It helps in
 !determining which J value is associated with which S1S2 term.
 !jval stores the total (2S1S2 type )coefficient associated with each J value
 !Each column represents a different J value for jval
 Real(kind = real_kind), Dimension(:,:,:), Allocatable :: jval_trkr
 !jval_trkr stores all the computed sets of J values
 Real(kind = real_kind) :: start, finish
 !start and finish help determine the time taken during the program execution
 Character (Len = 100) :: FileName1, FileName2, read_line, Outfile
 Logical :: EndOfFile, singular, normal, identical, stdv
 !Call cpu_time(start)
 !!Call init()

Contains
!!!Initialization Subroutines!!!
         Subroutine Init1()
           !!!!! Declaration and initailisation of local and global variables
           !!!!!!
                Logical :: check, check1
                Integer(kind = int_kind) :: pos
                !pos helps in defining the different J values.
                Real(kind = real_kind) :: val
                check = .False.
                check1 = .False.
                168 Format(i2)
                189 Format(F4.1)
                193 Format(F6.3)
                245 Format(F10.7)
                249 Format(F10.6)
           !!!!! Declaration and initailisation ends !!!!!

           !!!!! Determining number of magnetic centres (i.e. metal atoms) !!!!!
           !!!!! Initialisation of spin for each metal centre !!!!!
                Read(11, '(a)', iostat = ios) read_line
                If(trim(read_line) == 'magnetic centres') Then
                        Read(11, 168) num_mag_cent
                Else
                        Write(*,*) ' Erroneous number of magnetic centres in the input file'
                        Call Abort
                End If
                tot_interac = num_mag_cent * (num_mag_cent - 1) / 2
                Write(12,*) "Number of Magnetic centres:", num_mag_cent
                Write(12,*) "Total number of possible interactions:", tot_interac
           !!!!! Determination of magnetic centres and their initialization ends
           !!!!!!

           !!!!! Determining how many sets of spin densities are there for each
           !metal centre !!!!!
           !!!!! This is achieved by counting the number of lines in the spin
           !density file   !!!!!
                i = 0 !Basically will store the number of rows in the spin density file
                      !here.
                Write(12,*) "Input obtained from spin density file"
                Read(13, '(a)', iostat = ios) read_line
                Do While (EndofFile .ne. .True.)
                        If (ios .eq. Iostat_End) then
                                EndOfFile = .True.
                                Exit
                        End If
                        If (read_line(1:5) .ne. '     ' .AND. EndofFile .ne. .True.) Then
                                i = i + 1
                                Write (12,'(a)', advance = 'no') read_line
                                Write(12,*)
                        End If
                        Read(13, '(a)', iostat = ios) read_line
                End Do
                Rewind(13)
                Allocate (spin(i, num_mag_cent), energy(i))
                num_spin_dens_set = i
                Write (12, *) 'Number of spin density sets:', num_spin_dens_set
           !!!!! Determination Ends !!!!!

           !!!!! Reading in the actual values of spin for each centre !!!!!
           !!!!! Also reading in the energy for each configuration    !!!!!
                Do i = 1, num_spin_dens_set
                        Read (13,*) (spin(i,j),j=1, num_mag_cent), energy(i)
                End Do
           !!!!! Reading of spins and energies ends !!!!!

           !!!!! Printing the spin matrix !!!!!
                246 Format(2x, F10.7)
                116 Format(2x, F18.10)
                Write (12,*) 'Spin Matrix'
                Do i = 1, num_spin_dens_set
                        Do j = 1, num_mag_cent
                                spin(i,j) = nint(spin(i,j))
                                Write (12, 246, advance = 'no') spin(i,j)
                        End Do
                        Write (12, 116, advance = 'no') energy(i)
                        Write (12,*)
                End Do
           !!!!! Printing Ends !!!!!

           !!!!! Determining the number of J values, reading them and printing
           !them!!!!!
                Read(11, '(a)', iostat = ios) read_line
                If(trim(read_line) == 'J values') then
                        Read(11, 168) no_of_j_val
                Else
                        Write(*,*) ' Erroneous No. of J values'
                        Call Abort
                End If
                Allocate(jpos(no_of_j_val + 1))
                Allocate(stddev(Size(jpos) - 1), std_per(Size(jpos) - 1))
                Write (12,*)
                Write (12,*)'No. of J value asked for = ', no_of_j_val
                If (no_of_j_val >= num_spin_dens_set) Then
                        Write(12,*) ' You need more equations to obtain these many J values. '
                        Write(12,*) ' Aborting Run'
                        Call Abort()
                End If
                !The addition of 1 has been made to make keep the last value as
                !the total
                !number of interaction that will help in comparisons later
           !!!!! J value determination ends !!!!!

           !!!!! Reading in the Spin Hamiltonian !!!!!
                Read(11, '(a)') read_line
                If(trim(read_line) .ne. 'Hamiltonian') Then
                        Write(*,*) 'No Hamiltonian'
                        Call Abort
                End If
                Allocate (hamil1(tot_interac), hamil2(tot_interac), hamil(num_spin_dens_set,tot_interac))
                check = .False.
                pos = 2
                i = 1
                jpos(1) = 1
                hamil1 = 0
                hamil2 = 0
                Read(11, '(a)') read_line
                Do While (trim(read_line) .ne. 'Hamiltonian Ends')
                        If(trim(read_line) == '****') Then
                                jpos(pos) = i
                                pos = pos + 1
                        Else
                                Backspace(11)
                        End If
                        Read(11, *) hamil1(i), hamil2(i)
                        i = i + 1
                        Read(11, '(a)') read_line
                End Do
                jpos(pos) = i
                !The above line will store the total interaction in the last position of
                !jpos
                !Write(*,*) pos
                !If (pos < no_of_j_val) then
                !        Write(*,*)'You have specified more J values but
                !        actually have less'
                !        Write(*,*)'Aborting'
                !        Call abort()
                !End If

           !!!!! Calculating the weight of each interaction i.e., the term
           !-2S1S2
                Do i = 1, num_spin_dens_set
                        Do j = 1, tot_interac
                                hamil(i, j) = 0
                        End Do
                End Do
                Do i = 1, num_spin_dens_set
                        Do j = 1, tot_interac
                                If(hamil1(j) > 0 .And. hamil2(j) > 0) Then
                                        hamil(i,j) = -0.50 * spin(i, hamil1(j)) * spin(i, hamil2(j))
                !The above line is the same as 2S1S2 where Si = 0.5*(spin value
                !obtained from the density file)
                !Si has not been divided by 2 each time so the overall value has
                !been
                !multiplied by 2 in the end to take that part into account.
                                End If
                        End Do
                End Do
            !!!!! Printing the Hamiltonian !!!!!
                j = 2
                134 Format(i2)
                Write(12,*) 'Interactions considered under each J value'
                Write(12, '(a)', advance = 'no') 'J 1   '
                Do i = 1, tot_interac
                        If (hamil1(i) .ne. 0 .AND. hamil2(j) .ne. 0) then
                                If (i == jpos(j)) Then
                                        Write(12, *)
                                        j = j + 1
                                        Write(12, '(a)', advance = 'no') 'J'
                                        Write(12, 134, advance = 'no') (j-1)
                                        Write(12, '(a)', advance = 'no') '   '
                                End If
                        Write(12, 168, advance = 'no') hamil1(i)
                        Write(12, 168, advance = 'no') hamil2(i)
                        Write(12, '(a)', advance = 'no') '  '
                        End If
                End Do
                Write(12,*)
                !Debugger Write(12, *) 'jpos', jpos
                !Write(12, *) 'Coefficient associated with each interaction'
                !Write(12, '(a)', advance = 'no') ' '
                !Write(12, *) hamil1, hamil2
                !Do i = 1, num_spin_dens_set
                !        Do j = i+1, num_spin_dens_set
                !                If (hamil1(i) .ne. 0 .AND. hamil2(j) .ne. 0)
                !                then
                !                        Write(12, 168, advance='no') i
                !                        Write(12, 168, advance='no') j
                !                        Write(12, '(a)', advance = 'no') '
                !                        '
                !                End If
                !        End Do
                !End Do
                !Write(12, *)
                !Do i = 1, num_spin_dens_set
                !        Do j = 1, tot_interac
                !                If (hamil1(j) .ne. 0 .AND. hamil2(j) .ne. 0)
                !                then
                !                        Write(12, 249, advance='no') hamil(i,j)
                !                End If
                !        End Do
                !        Write(12, *)
                !End Do

         End Subroutine Init1



End Module init
