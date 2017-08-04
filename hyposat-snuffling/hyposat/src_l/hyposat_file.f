cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function file_check(name)
c
c     function checks if file exists and if eventually
c     the file must be searched under the directory
c     given by environment variable HYPOSAT_DATA
c
c     26 October 2000
c

      character*(*) name
      character*80 path, envvar, file, file_check
      integer trimle

      logical exst

      file_check = ''

      ifil = trimle(name)
      inquire(file=name(1:ifil),exist=exst)

      if(exst) then

         file_check = name(1:ifil)

      else

c
c        Start of a SUN-specific block
c        Get name of directory path from UNIX environment variable:
c        HYPOSAT_DATA. This part has eventually be changed for
c        other machines!!!
c
         envvar = 'HYPOSAT_DATA'
         CALL GETENV(ENVVAR, path)

         file = path(1:trimle(path)) // '/' // name(1:ifil)

         ifil = trimle(file)

         inquire(file=file(1:ifil),exist=exst)

         if(exst) then

            file_check = file(1:ifil)

         else

            write (*,'('' Major ERROR: Cannot find to open: '',a)') 
     *            name
            write (*,'('' Major ERROR: Cannot find to open: '',a)') 
     *            file

            write (*,'(''Set correct path name with environment '',
     *	               ''variable HYPOSAT_DATA !'')')
            stop

	 endif

      ENDIF

      return
      end
