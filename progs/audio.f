           implicit real*4 (a-h,o-z)

           character*6 HETATOM
           character*4 at
           character*3 RES
           dimension x(2000),y(2000),z(2000),nresat(2000)
           dimension contact(2000),contactacum(2000)
           logical hacer(2000)
           integer*4 final
            
           read(*,*) ta,d
           do i = 1,10000
           read(*,*) ta,d2
           d3 = (d-d2)*5000000*1
           final=int(d2*100)
            d=d2
           write(*,*) final, 0 
           enddo

           end

