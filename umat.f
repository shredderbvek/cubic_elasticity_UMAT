      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      integer i, j, k, l, m, n
      real C11, C12, C44
      real phi1, Phi, phi2
      real R(3,3), C(6,6), C_rot(6,6)
      real T(6,6), T_inv(6,6)
      real pi

      pi = 3.141592653589793

C -----------------------------------------------------------------
C Extract Material Properties
C -----------------------------------------------------------------
      C11 = props(1)
      C12 = props(2)
      C44 = props(3)

C Euler Angles in Degrees → Radians
      phi1 = props(4) * pi / 180.0
      Phi  = props(5) * pi / 180.0
      phi2 = props(6) * pi / 180.0

C -----------------------------------------------------------------
C Construct Stiffness Matrix for Cubic Symmetry (Voigt notation)
C -----------------------------------------------------------------
      do i = 1, 6
        do j = 1, 6
          C(i,j) = 0.0d0
        end do
      end do

      C(1,1) = C11
      C(2,2) = C11
      C(3,3) = C11

      C(1,2) = C12
      C(1,3) = C12
      C(2,1) = C12
      C(2,3) = C12
      C(3,1) = C12
      C(3,2) = C12

      C(4,4) = C44
      C(5,5) = C44
      C(6,6) = C44

C -----------------------------------------------------------------
C Compute Rotation Matrix (Bunge’s Convention)
C -----------------------------------------------------------------
      R(1,1) = COS(phi1)*COS(phi2) - SIN(phi1)*SIN(phi2)*COS(Phi)
      R(1,2) = -COS(phi1)*SIN(phi2) - SIN(phi1)*COS(phi2)*COS(Phi)
      R(1,3) = SIN(phi1)*SIN(Phi)

      R(2,1) = SIN(phi1)*COS(phi2) + COS(phi1)*SIN(phi2)*COS(Phi)
      R(2,2) = -SIN(phi1)*SIN(phi2) + COS(phi1)*COS(phi2)*COS(Phi)
      R(2,3) = -COS(phi1)*SIN(Phi)

      R(3,1) = SIN(phi2)*SIN(Phi)
      R(3,2) = COS(phi2)*SIN(Phi)
      R(3,3) = COS(Phi)

C -----------------------------------------------------------------
C Compute Transformation Matrix T (Voigt notation)
C -----------------------------------------------------------------
      call compute_voigt_transformation(R, T, T_inv)

C -----------------------------------------------------------------
C Rotate Stiffness Matrix
C -----------------------------------------------------------------
      do i = 1, 6
        do j = 1, 6
          C_rot(i,j) = 0.0d0
          do k = 1, 6
            do l = 1, 6
              C_rot(i,j) = C_rot(i,j) + T(i,k) * C(k,l) * T(j,l)
            end do
          end do
        end do
      end do

C -----------------------------------------------------------------
C Assign Rotated Stiffness Matrix to ddsdde
C -----------------------------------------------------------------
      do i = 1, ntens
        do j = 1, ntens
          ddsdde(i,j) = C_rot(i,j)
        end do
      end do

C -----------------------------------------------------------------
C Stress Update Using Rotated Stiffness Matrix
C -----------------------------------------------------------------
      do i = 1, ntens
        do j = 1, ntens
          stress(i) = stress(i) + ddsdde(i,j) * dstran(j)
        end do
      end do

      return
      end


C -----------------------------------------------------------------
C Subroutine: Compute Voigt Transformation Matrix
C -----------------------------------------------------------------
      subroutine compute_voigt_transformation(R, T, T_inv)
      real R(3,3), T(6,6), T_inv(6,6)
      integer i, j

C Initialize matrices
      do i = 1, 6
        do j = 1, 6
          T(i,j) = 0.0d0
          T_inv(i,j) = 0.0d0
        end do
      end do


      T(1,1) = R(1,1)**2
      T(1,2) = R(1,2)**2
      T(1,3) = R(1,3)**2
      T(1,4) = 2*R(1,2)*R(1,3)
      T(1,5) = 2*R(1,3)*R(1,1)
      T(1,6) = 2*R(1,1)*R(1,2)
      T(2,1) = R(2,1)**2
      T(2,2) = R(2,2)**2
      T(2,3) = R(2,3)**2
      T(2,4) = 2*R(2,2)*R(2,3)
      T(2,5) = 2*R(2,3)*R(2,1)
      T(2,6) = 2*R(2,1)*R(2,2)
      T(3,1) = R(3,1)**2
      T(3,2) = R(3,2)**2
      T(3,3) = R(3,3)**2
      T(3,4) = 2*R(3,2)*R(3,3)
      T(3,5) = 2*R(3,3)*R(3,1)
      T(3,6) = 2*R(3,1)*R(3,2)

      T(4,1) = R(2,1)*R(3,1)
      T(4,2) = R(2,2)*R(3,2)
      T(4,3) = R(2,3)*R(3,3)
      T(4,4) = R(2,2)*R(3,3) + R(2,3)*R(3,2)
      T(4,5) = R(2,1)*R(3,3) + R(2,3)*R(3,1)
      T(4,6) = R(2,2)*R(3,1) + R(2,1)*R(3,2)
      T(5,1) = R(1,1)*R(3,1)
      T(5,2) = R(1,2)*R(3,2)
      T(5,3) = R(1,3)*R(3,3)
      T(5,4) = R(1,2)*R(3,3) + R(1,3)*R(3,2)
      T(5,5) = R(1,3)*R(3,1) + R(1,1)*R(3,3)
      T(5,6) = R(1,1)*R(3,2) + R(1,2)*R(3,1)
      T(6,1) = R(1,1)*R(2,1)
      T(6,2) = R(1,2)*R(2,2)
      T(6,3) = R(1,3)*R(2,3)
      T(6,4) = R(1,2)*R(2,3) + R(1,3)*R(2,2)
      T(6,5) = R(1,3)*R(2,1) + R(1,1)*R(2,3)
      T(6,6) = R(2,2)*R(1,1) + R(2,1)*R(1, 2)

      return
      end
