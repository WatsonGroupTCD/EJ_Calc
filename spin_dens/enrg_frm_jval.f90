 !The subroutines that can be called from here are avg_spins and all_poss_spins
 !All the others are subroutines to support these two subroutines

Module enrg_frm_jval 
 Use Iso_Fortran_Env
 Use init
 Implicit None
 Real(kind = real_kind), Dimension(:), Allocatable :: unp_elec
 !unp_elec essentially stores the number of enpaired electrons on each metal
 !centre
 Real(kind = real_kind), Dimension(:,:), Allocatable :: all_config
 !all_config stores all the possible spin configurations of spin
 Real(kind = real_kind), Dimension(:), Allocatable:: new_nrg
 !new_nrg stores the energy of all possible spin configurations
 Integer(kind = int_kind) :: poss_comb1
 !poss_comb1 stores the total possible spin configurations

Contains
!The following subroutine stores the average spin on each metal centre. This is
!calculated by taking the absolute average of all formal spin values that are
!known for the different configurations (in the spin density file).  
Subroutine avg_spin()
 Integer (kind = int_kind) :: i1, j1
 poss_comb1 = 2**(num_mag_cent-1)
 Allocate (unp_elec(num_mag_cent), all_config(num_mag_cent, poss_comb1))
 Allocate (new_nrg(poss_comb1))
 unp_elec = 0
 !Determining the average spin value on each metal centre in the ferromagnetic
 !configuration
 Do i1 = 1, num_mag_cent
        Do j1 = 1, num_spin_dens_set
                unp_elec(i1) = unp_elec(i1) + Abs(spin(j1,i1))
        End Do
 End Do
 unp_elec = unp_elec/num_spin_dens_set
 Write(16,*) 'This file contains the energy (in cm-1)  of all the possible spin'
 Write(16,*) 'configurations using the J value set obtained.'
 Write(16,*) ' The J value set (cm-1) used is as follows:'
 Write(16,*) tot_avg
 Write(16,*) 'Total possible spin configurations:', poss_comb1
 Write(16,*) 'Average Spin Density on each metal centre:'
 118 Format (1x, F8.5) 
 Do i1 = 1, num_mag_cent
        Write(16, 118, advance = 'no') unp_elec(i1)
 End Do 
 Call std_dev2d(unp_elec, spin)
End Subroutine avg_spin

!The following subroutine alongwith the subroutine update_spin forms all
!possible spin configurations, stores them in all_config, determines the energy
!of each configuration in accordance with the specified hamiltonian  and stores
!the energy in new_nrg
Subroutine all_poss_spins()
 Integer (kind = int_kind) :: tot_atoms, i1, j1, k1
 115 Format (i4)
 119 Format (1x, F12.4)
 118 Format (1x, F8.5)
 116 Format (i3)
 Do j1 = 1, num_mag_cent 
        Write(16, '(a)', advance = 'no') ' M'
        Write(16, 115, advance = 'no') j1
        Write(16, '(a)', advance = 'no') '      '
 End Do
 Write(16, '(a)', advance = 'no') '    Energy(cm-1)'
 Write(16, *)
 new_nrg = 0
 all_config = 0
 Do i1 = 1, poss_comb1 
        k1 = 1
        all_config(:,i1) = unp_elec
        Do j1 = 1, Size(hamil1)
               If (j1 == jpos(k1 + 1)) then
                       k1 = k1 + 1
               End If
               new_nrg(i1) = new_nrg(i1) - 0.5*tot_avg(k1)*unp_elec(hamil1(j1))*unp_elec(hamil2(j1))
        End Do
        Do k1 = 1, num_mag_cent
                Write(16, 118, advance = 'no') all_config(k1,i1)
                Write(16, '(a)', advance = 'no') '   '
        End Do
        Write(16, 119) new_nrg(i1)
        If (unp_elec(num_mag_cent) > 0) then
                unp_elec(num_mag_cent) = unp_elec(num_mag_cent) * -1
        Else
                Call update_spin()
        End If
 End Do
 Write(16, *) 'Sorted Energies (Energies (cm-1) wrt to the ground state):'
 Do j1 = 1, num_mag_cent
        Write(16, '(a)', advance = 'no') ' M'
        Write(16, 115, advance = 'no') j1
        Write(16, '(a)', advance = 'no') '      '
 End Do
 Write(16, '(a)', advance = 'no') '    Energy(cm-1)'
 Write(16, *)
 Call sort()
End Subroutine all_poss_spins

!Bubble sorting
Subroutine sort()
 Integer (kind = int_kind) :: tot_atoms, i1, j1, k1
 Do i1 = 1, poss_comb1
        Do j1 = 1, poss_comb1 - i1 
                If (new_nrg(j1) > new_nrg(j1 + 1)) then
                        Call swap(j1)
                End If
        End Do
 End Do
 Call prnt() 
End Subroutine sort

!Helper routine for the swap subroutine. It swaps two different values
Subroutine swap(a)
 Integer (kind = int_kind) :: a
 Real (kind = real_kind) :: nrg_swp
 Real (kind = real_kind), Dimension(num_mag_cent) :: unp_elec_swp
 nrg_swp = new_nrg(a)
 new_nrg(a) = new_nrg(a + 1)
 new_nrg(a + 1) = nrg_swp
 unp_elec_swp = all_config(:,a)
 all_config(:,a) = all_config(:, a+ 1)
 all_config(:, a + 1) = unp_elec_swp
 
End Subroutine swap

!Subroutine to print the sorted energies in the nrg_frm_jval file
Subroutine prnt()
 Integer (kind = int_kind) :: i1, k1
 119 Format (1x, F12.4)
 116 Format (i3)
 118 Format (1x, F8.5)
 Do i1 = 1, poss_comb1 - 1
        Do k1 = 1, num_mag_cent
                Write(16, 118, advance = 'no') all_config(k1,i1)
                Write(16, '(a)', advance = 'no') '   '
        End Do
        Write(16, 119) new_nrg(i1) - new_nrg(1)
        If (Abs(new_nrg(i1) - new_nrg(i1 + 1)) > 1.0) Then
                Write(16,*)
        End If
 End Do
 Do k1 = 1, num_mag_cent
        Write(16, 118, advance = 'no') all_config(k1,i1)
        Write(16, '(a)', advance = 'no') '   '
 End Do
 Write(16, 119) new_nrg(i1)  - new_nrg(1)
 Write(16,*) ' ---End Of File---'
End Subroutine prnt

!helper subroutine for the subroutine all_poss_spins. This subroutine allows
!going through all possible spin configurations without going out of bounds.
Subroutine update_spin()
 Integer (kind = int_kind) :: i2, j2
 Do i2 = num_mag_cent, 1, -1
         If (unp_elec(i2) > 0) then
                 unp_elec(i2) = unp_elec(i2) * -1
                 Do j2 = i2 + 1, num_mag_cent
                         unp_elec(j2) = Abs(unp_elec(j2))
                 End Do
                 exit
         End If
 End Do
        
End Subroutine update_spin

Subroutine std_dev2d(avg, full_dat)
       Integer(kind = int_kind) :: k, l, m
       Real(kind = real_kind):: ss
       Real(kind = real_kind), Dimension(:), Allocatable :: avg, var, stddev, std_per
       Real(kind = real_kind), Dimension(:,:), Allocatable :: full_dat
       Allocate (var(Size(avg)), stddev(Size(avg)), std_per(Size(avg)))
       var = 0
       stddev = 0
       Do i = 1, Size(full_dat)/Size(avg)
               Do j = 1, Size(avg)
                       var(j) = var(j) + (avg(j) - Abs(full_dat(i, j)))**2
               End Do
       End Do
       var = var/main_cntr
       stddev(:) = sqrt(var(:))
       std_per(:) = stddev(:) / Abs(avg(:)) * 100
       Write(16,*)'standard deviation for spins in the spin density file:'
       165 Format(1x, F7.4)
       166 Format(i3)
       Do j = 1, Size(avg)
               Write(16, 165, advance = 'no') stddev(j)
               Write(16, '(a)', advance = 'no') '('
               Write(16, 166, advance = 'no') Int(std_per(j))
               Write(16, '(a)', advance = 'no') '%)'
       End Do
       Write(16,*)
       !Call std_dev1(stddev)
End Subroutine std_dev2d

End Module enrg_frm_jval

