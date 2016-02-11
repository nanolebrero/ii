
           implicit real*8 (a-h,o-z)

            real*8 , dimension (:), allocatable :: x,y
            integer*4 listo
            real*8, dimension (:), allocatable :: fx,fy,vx,vy
            open(unit=56,file='geom.xyz')
            open(unit=23,file='salida')
            open(unit=55,file='ini')
           read(55,*) npas,npart,fact,fact2,dxmax
           allocate(x(npart),y(npart),fx(npart),fy(npart)
     > ,vx(npart),vy(npart))
           write(*,*) npart
c           do jj=1,npart
c           read(55,*) x(jj),y(jj)
c           enddo
            kwrite=0            
            x=0
            vy=0
            vx=0
           ! calculo fuerzas
            do j=1,npas
            fx=0
c            fy=0
           do i=1,npart-1
              dx=x(i+1)-x(i)

c              dx=min(dx,dxmax)
c              dx2=dx*dx
              fx(i)=fx(i) + dx !+ 0.0005*dx2
              fx(i+1)=fx(i+1) -dx !- 0.0005*dx2
           enddo
           if((j.gt.100.and.j.lt.400).or.(j.gt.200011.and.j.lt.200325))
     >      then
           do kk= 1,npart
            fx(kk)=fx(kk)+5.D-3*exp(-0.0001*((kk-3*npart/10)**2))
           enddo

           endif
            
           do i=1,100
            fx(i)=fx(i) - x(i)*exp(-0.0001*real(i))*10.1
c            write(*,*) 'aca estoy'
            fx(npart-i)=fx(npart-i)-x(npart-i)*exp(-0.0001*real(10-i))*10.1
           enddo           


           do i=1,npart
            vx(i)=vx(i)+fx(i)*fact2
            vx(i) = vx(i)*fact
           enddo

            do i=1,npart
c            write(44,*) j, i, vx(i) 
            x(i)=x(i)+vx(i)
            enddo
             x(1)=0.
             x(npart)=0.
             if (kwrite.eq.20) then
              listo=100000*x(npart/4)
               kwrite=0
               write(23,*) listo 
             else
              kwrite=kwrite+1
            endif
            
           enddo
           end

