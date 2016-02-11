
           implicit real*4 (a-h,o-z)

            real*4 , dimension (:), allocatable :: x,y
            integer*4 listo
            real*4, dimension (:), allocatable :: fx,fy,vx,vy
            logical toco 
            open(unit=56,file='geom.xyz')
            open(unit=23,file='salida')
            open(unit=55,file='ini')
           read(55,*) npas,npart,fact0,fact2,dxmax
           allocate(x(npart),y(npart),fx(npart),fy(npart)
     > ,vx(npart),vy(npart))
           write(*,*) npart
            kwrite=0            
            x=0
            vy=0
            vx=0
           ! calculo fuerzas
              
            do i=1,npart
            fy(i) = 1.D-5*exp(-0.001*
     >      ((i-npart/4)**2))
            enddo


            do j=1,npas
            fx=0
c            fy=0
            do i=2,npart-1
              fx(i)=x(i-1)+x(i+1)-2*x(i) !+ 0.0005*dx2
            enddo
c            fx(npart)=x(npart)
c            if(j.eq.2000000) fact2=fact2*2
           toco=.false.
           if((j.gt.100.and.j.lt.400000).
     >     or.(j.gt.800000.and.j.lt.1500000). 
     >     or.(j.gt.1600000.and.j.lt.2500000). 
     >     or.(j.gt.2900000.and.j.lt.3500000)) 
     >      then
           toco=.true. 
          if(mod(j,2500).eq.0) then
            azar=rand()-0.5

           fff=real(j)/1000.

           fff=min(1.,fff)
           azar = fff*azar
           endif

          if(mod(j,2500).lt.250) then
c           do kk= 1,npart
           fx(npart-1)=fx(npart-1)+azar
c           enddo

           endif
           endif
           if (toco) then
           fact=fact0
           else
            fact=fact0*100
            endif
           do i=1,npart
            vx(i)=vx(i)+fx(i)*fact2
            vx(i) = vx(i)*(1.-fact)
           enddo

            do i=1,npart
c            write(44,*) j, i, vx(i) 
            x(i)=x(i)+vx(i)
            enddo
             x(1)=0.
             x(npart)=0.
             listo=0
             alisto=0
             if (kwrite.eq.1) then

        
              do ii=1,25
              alisto=alisto + x(ii)
               enddo
c              alisto=x(1)-x(2)
              listo=alisto*10000
               kwrite=0
               write(23,*) listo 
             else
              kwrite=kwrite+1
            endif
            
           enddo
           end

