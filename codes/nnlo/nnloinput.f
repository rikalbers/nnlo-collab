      function nnloinput(stringa)
      use process
      implicit none
      real * 8 nnloinput
      character (len=*) stringa
      integer maxnum
      parameter (maxnum=100)
      character * 100 fname
      character * 100 line,line0
      character * 20 string
      character * 20 keywords(maxnum)
      real * 8 values(maxnum)
      integer ios,numvalues,j,k,l,imode
      integer ini
      data ini/0/      
      save ini, keywords, values, numvalues, fname
      string = stringa
      if(ini.eq.0) then
         fname = process_name(1:len_trim(process_name))//'-nnlo.input'
         open(unit=33,file=fname,status='old',iostat=ios)
         if (ios.ne.0) then
           print *,"Problem occured with input card..."
           print *,fname(1:len_trim(fname))," is not present..."
           stop
         end if
         numvalues=0
         do l=1,1000000
            line0=' '
            read(unit=33,fmt='(a)',iostat=ios) line0
            if(ios.ne.0.and.line0.eq.' ') goto 10
            line=line0
            do k=1,100
               if(line(k:k).eq.'#'.or.line(k:k).eq.'!') then
                  line(k:)=' '
               endif
            enddo
            if(line.ne.' ') then
               if(numvalues.eq.maxnum) then
                  write(*,*) ' too many entries in nnloinput.dat'
                  stop
               endif
               numvalues=numvalues+1
c skip blanks
 12            if(line(1:1).eq.' ') then
                  line=line(2:)
                  goto 12
               endif
               k=index(line,' ')
               keywords(numvalues)=line(1:k-1)
               line=line(k+1:)
               read(unit=line,fmt=*,iostat=ios) values(numvalues)
               if(ios.ne.0) then
                  write(*,*) ' nnloinput error: cannot parse '
                  write(*,'(a)') line0
                  stop
               endif
            endif
         enddo
 10      continue
         close(33)
         ini=1
      endif
      if(string(1:1).eq.'#') then
         string=string(2:)
         imode=0
      else
         imode=1
      endif
      do j=1,numvalues
         if(string.eq.keywords(j)) then
            nnloinput=values(j)
            return
         endif
      enddo
      if(imode.eq.1) then
         write(*,*) ' nnloinput: keyword ',string,' not found'
         stop
      endif
      nnloinput=-1d6
      end
