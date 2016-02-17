
           implicit real*4 (a-h,o-z)

            real*4 , dimension (:,:), allocatable :: x,y
            integer*4 listo
            real*4, dimension (:,:), allocatable :: fx,fy,vx,vy
            logical toco 
            open(unit=56,file='geom.xyz')
            open(unit=23,file='salida')
            open(unit=55,file='ini')
           read(*,*) npas,npart,fact0,fact2,dxmax
           allocate(x(npart,npart),fx(npart,npart),fy(npart,npart)
     > ,vx(npart,npart),vy(npart,npart))
           write(*,*) npart
            kwrite=0            
            x=0
            vy=0
            vx=0

              !!!!!!!!!!!! Precalculo de gaussiana de exitacion
            do i=1,npart
             do ii=1,npart
            fy(i,ii) = 2.D-6*exp(-0.0001*
     >      ((i-npart/4)**2))*exp(-0.0001*((ii-npart/4)**2))
            enddo
            enddo
              !!!!!!!!!!!!!!!
            do j=1,npas
            fx=0
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculo fuerzas
            do ii=2,npart-1
            do i=2,npart-1
             fx(i,ii)=x(i-1,ii)+x(i+1,ii)-4*x(i,ii)+x(i,ii-1)+x(i,ii+1) !+ 0.0005*dx2
              
c              if(fx(i).ge.0) then
c              fx(i)= min(fx(i),dxmax)
c              else
c              fx(i)=max(fx(i),-1*dxmax)
c              endif
             enddo
            enddo
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !!!!!!!!!!!!exitación
           toco=.true.
           if((j.gt.100.and.j.lt.400000).
     >     or.(j.gt.800000.and.j.lt.1500000). 
     >     or.(j.gt.1600000.and.j.lt.2500000). 
     >     or.(j.gt.2900000.and.j.lt.3500000)) 
     >      then
           toco=.true. 
          if(mod(j,2500).eq.0) then !!!El valor aleatorio lo actualizo cada tanto para evitar ruido muy agudo
            azar=rand()-0.5
           fff=real(j)/1000.
           fff=min(1.,fff)
           azar = fff*azar
           endif

          if(mod(j,250000).lt.250) then   !!!! El valor aleatorio lo aplico 1/10 del tiempo
           do kkk= 1,npart
            do kk=1,npart
           fx(kk,kkk)=fx(kk,kkk)+fy(kk,kkk)!*azar    !!!! Aplicando
           enddo
           enddo
           endif
           endif
           if (toco) then   
           fact=fact0         !!!!! Si estoy tocando el pagamiento disminuye
           else
            fact=fact0*100    !!!!! Si dejo de tocar el apagamiento aumenta
            endif
            do ii=1,npart
           do i=1,npart
            vx(i,ii)=vx(i,ii)+fx(i,ii)*fact2  !!!! Actualizao velocidades por la fuerza
            vx(i,ii) = vx(i,ii)*(1.-fact)  !!!! Aplico fricción 
            enddo
           enddo
            npartm=npart/2
            npartm2=(npartm-1)**2 

 
             do ii=1,npart
            do i=1,npart
c            write(44,*) j, i, vx(i) 
            x(i,ii)=x(i,ii)+vx(i,ii) !!!! Actualizo posiciones
             iquieto=(i-npartm)**2+(ii-npartm)**2
            if(iquieto.ge.npartm2) then
              x(i,ii)=0
c             write(*,*) i, ii

            endif
c            if(i.ge.npart.or.ii.ge.npart) x(i,ii)=0.
              

            enddo
           enddo
c              stop

c             x(1)=0.    !!!! Los extremos estan quietos
c             x(npart)=0.

             listo=0   !!!!Output en integer*4
             alisto=0  !!!! output en Real*4
             if (kwrite.eq.2) then  !!! Cada cuanto escribe
              do ii=npart/4-5,npart/4+5 !!! Acumulo sobre varios puentos para filtrar ultrasonido
                do i=npart/4-5,npart/4+5
              alisto=alisto + x(ii,i)
                enddo
               enddo
              listo=alisto*1000000 !!!! Multiplico por un factor que de numeros razonables para el kwave
               kwrite=0
               write(23,*) listo 
             else
              kwrite=kwrite+1
            endif
            
           enddo
           end

