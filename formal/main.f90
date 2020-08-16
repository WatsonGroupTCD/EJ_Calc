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



!This program calculates the J values using the formal spin density irrespective
!of the type of spin density provided in the spin density file. It converts the
!spin density to integer values by rounding them to the closest integer.
!Additionally, this program removes all singular solutions i.e. the solutions in
!which at least one of the equations becomes completely zero.


!It is recommended that you run this program in conjunction with the next
!version because that uses the spin density values which are provided in the
!spin density file but at the same time does not calculate the set of equations
!that were found to be singular by this code.

!All subroutines called here are in the modul_lib.f90 file unless specified
!otherwise

Program jval_calc_form_comp
 Use Iso_Fortran_Env
 Use init
 Use modul_lib 
 Use enrg_frm_jval
 Call Get_Command_Argument(1, FileName1)
 Call Get_Command_Argument(2, FileName2)
 Outfile = trim(Filename1) // '_form.txt'
 stdv = .FALSE.
 Open (Unit = 11, file = FileName1, action = 'read', position = 'rewind', iostat = ios)
 Open (Unit = 12, file = Outfile, status = 'unknown', action = 'readwrite')
 Open (Unit = 13, file = FileName2, action = 'read', position = 'rewind', iostat = ios)
 Open (Unit = 14, file = 'singul', status = 'unknown', action = 'readwrite')
 Open (Unit = 15, file = 'bigg', status = 'unknown', action = 'readwrite')
 Open (Unit = 16, file = 'nrg_frm_jval', status = 'unknown', action = 'readwrite')
 Call cpu_time(start)
 Write(*,*) 'Report wiil be available in the file', Outfile
 Write(*,*) 'Energy related to energy of all possible spin configs'
 Write(*,*) 'is available in the file nrg_frm_jval'
 Write(12,*) 'Please cite: Altering the nature of coupling by changing the oxidation state in a {Mn6} cage, DOI: 10.1039/D0DT01404D.'
 Write(12,*)
 Write(12,*) '-------Calculating J values from Formal Spins-------'
 Write(12,*) 'This program will convert spin denstities to formal spins no'
 Write(12,*) 'matter what you give and will use the formal spins to solve everything.'
 Write(12,*) 'If you want to use spin densities, use a different code. Two extra'
 Write(12,*) "files 'bigg' and 'singul' will also be written. 'singul' is required"
 Write(12,*) "for any code the code that uses the spin densities. The file 'bigg'"
 Write(12,*) " tells you which set of equation has large values (>50 cm-1)"
 Write(12,*)
 
 Call init1() !from init.f90
 Call jvals() !Setting up the J value matrix (jval)  and determining the coefficients of
              !each J value
 Call solv() !Obtaining all possible J values.
 Call avrg_calcs()
 !The following loop section puts a constraint on the J value sets that will be
 !considered by choosing only those sets that are within 3 standard deviations.
 !The loop essentially checks the standar deviation then removes the values
 !beyond 3 standard deviations and recalculate the standard deviation and repeat
 !this cycle till a self consistency is achieved.
 139 Format(i2)
 Write(12,*) 'Removing solutions in which J-values deviate by more than 3 standard deviations'
 Write(12,*)'Standard deviation (percentage standard deviation) on different J-values:'
 Write (12,'(a)', advance = 'no')  'non-singular equations'
 Do i = 1, Size(jpos) - 1
         Write (12,'(a)', advance = 'no') '      '
         Write (12,'(a)', advance = 'no') 'J'
         Write (12, 139, advance = 'no') i
         Write (12,'(a)', advance = 'no') '    '
 End Do
 
 Do While(stdv == .FALSE.)
        Call std_dev(tot_avg)
 End Do
 Write(12,*)
 Write(12,*) 'Self-consistency achieved!'
 Call final_print()
 Call backtrack1(jval, tot_avg)
 Call avg_spin()       !from enrg_frm_jval.f90
 Call all_poss_spins() !from enrg_frm_jval.f90
 Call cpu_time(finish)
 Write(12,*) 'Time taken(sec) in program execution', finish-start
 !Write(12,*) 'Report available in the file', Outfile !Does not work for some reason
 Write(12,*) '---End Of File---'
 Deallocate (jpos, hamil, hamil1, hamil2, spin, jval_trkr)
 Close(11)
 Close(12)
 Close(13)
 Close(14)
 Close(15)
 Close(16)
End Program jval_calc_form_comp
