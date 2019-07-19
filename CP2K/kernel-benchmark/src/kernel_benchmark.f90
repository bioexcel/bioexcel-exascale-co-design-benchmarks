PROGRAM kernel_benchmark
  USE cell_types,                      ONLY: cell_type,&
                                             cell_create
  USE kinds,                           ONLY: dp
  USE qmmm_gpw_forces,                 ONLY: qmmm_forces_with_gaussian_LG
  
  IMPLICIT NONE

  ! variables that will be passed to kernel
  INTEGER                                            :: pgfs_size  
  REAL(KIND=dp), DIMENSION(3)                        :: cgrid_pw_grid_dr
  REAL(KIND=dp)                                      :: cgrid_pw_grid_dvol
  INTEGER, DIMENSION(2, 3)                           :: cgrid_pw_grid_bounds
  INTEGER, DIMENSION(2, 3)                           :: cgrid_pw_grid_bounds_local
  REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cgrid_pw_grid_cr3d
  INTEGER                                            :: num_mm_atoms
  REAL(KIND=dp), DIMENSION(:), POINTER               :: mm_charges 
  INTEGER, DIMENSION(:), POINTER                     :: mm_atom_index
  REAL(KIND=dp), DIMENSION(:,:), POINTER             :: mm_particles_r  
  INTEGER                                            :: para_env_num_pe
  INTEGER                                            :: para_env_mepos
  REAL(KIND=dp), DIMENSION(:, :), POINTER            :: Forces
  REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: per_pot_cr3d
  INTEGER, DIMENSION(3)                              :: per_pot_npts
  REAL(KIND=dp), DIMENSION(3)                        :: per_pot_dr
  INTEGER, DIMENSION(:), POINTER                     :: per_pot_mm_atom_index
  TYPE(cell_type), POINTER                           :: mm_cell
  REAL(KIND=dp), DIMENSION(3)                        :: dOmmOqm
  INTEGER                                            :: iw, par_scheme
  REAL(KIND=dp), DIMENSION(2)                        :: qmmm_spherical_cutoff
  LOGICAL                                            :: shells
  
  ! variables used in wrapper only
  INTEGER                                            :: mm_cell_id_nr, mm_cell_ref_count, mm_cell_symmetry_id
  LOGICAL                                            :: mm_cell_orthorhombic
  REAL(KIND=dp)                                      :: mm_cell_deth
  INTEGER, DIMENSION(3)                              :: mm_cell_perd
  REAL(KIND=dp), DIMENSION(3, 3)                     :: mm_cell_hmat, mm_cell_h_inv
  INTEGER                                            :: myunit
  INTEGER                                            :: i,j,k
  INTEGER                                            :: imax, jmax, kmax
  
  
  ! read pgfs_size
  OPEN(newunit=myunit, file='./data/pgfs_size.dat', action='READ', status='OLD')
  READ(myunit, *) pgfs_size
  CLOSE(myunit)


  ! read cgrid
  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%dr.dat', action='READ', status='OLD')
  READ(myunit, *) cgrid_pw_grid_dr
  CLOSE(myunit)
  
  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%dvol.dat', action='READ', status='OLD')
  READ(myunit, *) cgrid_pw_grid_dvol
  CLOSE(myunit)

  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%bounds.dat', action='READ', status='OLD')
  DO i=1,2
     DO j=1,3
        READ(myunit, *) cgrid_pw_grid_bounds(i,j)
     END DO
  END DO
  CLOSE(myunit)
  
  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%bounds_local.dat', action='READ', status='OLD')
  DO i=1,2
     DO j=1,3
        READ(myunit, *) cgrid_pw_grid_bounds_local(i,j)
     END DO
  END DO
  CLOSE(myunit)

  OPEN(newunit=myunit, file='./data/cgrid%pw_grid%dr.dat', action='READ', status='OLD')
  READ(myunit, *) cgrid_pw_grid_dr
  CLOSE(myunit)
  
  OPEN(newunit=myunit, file='./data/cgrid%cr3d.dat', action='READ', status='OLD')
  READ(myunit, *) imax
  READ(myunit, *) jmax
  READ(myunit, *) kmax
  ALLOCATE (cgrid_pw_grid_cr3d(imax,jmax,kmax))
  DO k=1,kmax
     DO j=1,jmax
        DO i=1,imax
           READ(myunit, *) cgrid_pw_grid_cr3d (i,j,k)
        END DO
     END DO
  END DO
  CLOSE(myunit)


  ! read num_mm_atoms
  OPEN(newunit=myunit, file='./data/num_mm_atoms.dat', action='READ', status='OLD')
  READ(myunit, *) num_mm_atoms
  CLOSE(myunit)

  ! read mm_charges
  OPEN(newunit=myunit, file='./data/mm_charges.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (mm_charges(imax))
  READ(myunit, *) mm_charges
  CLOSE(myunit)

  ! read mm_atom_index
  OPEN(newunit=myunit, file='./data/mm_atom_index.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (mm_atom_index(imax))
  READ(myunit, *) mm_atom_index
  CLOSE(myunit)
  
  ! read mm_particles
  OPEN(newunit=myunit, file='./data/mm_particles%r.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (mm_particles_r(imax,3))
  DO i=1,SIZE(mm_particles_r)
     READ(myunit, *) mm_particles_r(i,:)
  END DO
  CLOSE(myunit)
  
  ! read para_env
  OPEN(newunit=myunit, file='./data/para_env.dat', status='OLD', action='READ')
  WRITE(myunit, *) para_env_num_pe
  WRITE(myunit, *) para_env_mepos
  CLOSE(myunit)
  
  ! read Forces
  OPEN(newunit=myunit, file='./data/Forces.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  READ(myunit, *) jmax
  ALLOCATE (Forces(imax,jmax))
  DO j=1,jmax
     DO i=1,imax
        READ(myunit, *) Forces
     END DO
  END DO
  CLOSE(myunit)


  ! read per_pot
  OPEN(newunit=myunit, file='./data/per_pot%TabLR%cr3d.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  READ(myunit, *) jmax
  READ(myunit, *) kmax
  ALLOCATE (per_pot_cr3d(imax,jmax,kmax))
  DO k=1,kmax
     DO j=1,jmax
        DO i=1,imax
           READ(myunit, *) per_pot_cr3d(i,j,k) 
        END DO
     END DO
  END DO
  CLOSE(myunit)

  OPEN(newunit=myunit, file='./data/per_pot%TabLR%pw_grid%npts.dat', status='OLD', action='READ')
  READ(myunit, *) per_pot_npts
  CLOSE(myunit)
  
  OPEN(newunit=myunit, file='./data/per_pot%TabLR%pw_grid%dr.dat', status='OLD', action='READ')
  READ(myunit, *) per_pot_dr
  CLOSE(myunit)

  OPEN(newunit=myunit, file='./data/per_pot%mm_atom_index.dat', status='OLD', action='READ')
  READ(myunit, *) imax
  ALLOCATE (per_pot_mm_atom_index(imax))
  READ(myunit, *) per_pot_mm_atom_index
  CLOSE(myunit)
  
  ! read mm_cell
  OPEN(newunit=myunit, file='./data/mm_cell.dat', status='OLD', action='READ')
  READ(myunit, *) mm_cell_id_nr
  WRITE(myunit, *) mm_cell_ref_count
  WRITE(myunit, *) mm_cell_symmetry_id
  WRITE(myunit, *) mm_cell_orthorhombic
  WRITE(myunit, *) mm_cell_deth
  WRITE(myunit, *) mm_cell_perd
  WRITE(myunit, *) mm_cell_hmat
  WRITE(myunit, *) mm_cell_h_inv
  CLOSE(myunit)
  
  CALL cell_create(mm_cell,&
       mm_cell_id_nr, mm_cell_ref_count, mm_cell_symmetry_id,&
       mm_cell_orthorhombic,&
       mm_cell_deth, mm_cell_hmat, mm_cell_h_inv, mm_cell_perd)
  
    
  ! read d0mm0qm
  OPEN(newunit=myunit, file='./data/dOmmOqm.dat', status='OLD', action='READ')
  READ(myunit, *) dOmmOqm
  CLOSE(myunit)
  
  
  ! read iw
  OPEN(newunit=myunit, file='./data/iw.dat', status='OLD', action='READ')
  READ(myunit, *) iw
  CLOSE(myunit)
  
  
  ! read par_scheme
  OPEN(newunit=myunit, file='./data/par_scheme.dat', status='OLD', action='READ')
  READ(myunit, *) par_scheme
  CLOSE(myunit)
  
  
  ! read qmmm_spherical_cutoff
  OPEN(newunit=myunit, file='./data/qmmm_spherical_cutoff.dat', status='OLD', action='READ')
  READ(myunit, *) qmmm_spherical_cutoff
  CLOSE(myunit)
  
  
  ! read shells
  OPEN(newunit=myunit, file='./data/shells.dat', status='OLD', action='READ')
  READ(myunit, *) shells
  CLOSE(myunit)
      
  
       
      
      
      
  
      
  
      
  ! CALL qmmm_forces_with_gaussian_LG  (pgfs_size,&
  !       cgrid_pw_grid_dr, cgrid_pw_grid_dvol, cgrid_pw_grid_bounds, cgrid_pw_grid_bounds_local, cgrid_pw_grid_cr3d,&
  !       num_mm_atoms,&
  !       mm_charges,& 
  !       mm_atom_index,&
  !       mm_particles_r,&
  !       para_env_num_pe, para_env_mepos,&
  !       Forces,&
  !       per_pot_cr3d, per_pot_npts, per_pot_dr, per_pot_mm_atom_index,&
  !       mm_cell,&
  !       dOmmOqm,&
  !       iw,&
  !       par_scheme,&
  !       qmmm_spherical_cutoff,&
  !       shells)



END PROGRAM
