      SUBROUTINE TANG(AS,BS,NUUR)
      IMPLICIT NONE
#include "truncation.h"
      REAL*8 AS(NSH2,3),BS(NSH2,3)
      INTEGER NUUR
      print *,"FAKE TANGENT LINEAR"
      BS=AS
      RETURN
      END

      SUBROUTINE ADJ(AS,BS,NUUR)
      IMPLICIT NONE
#include "truncation.h"
      REAL*8 AS(NSH2,3),BS(NSH2,3)
      INTEGER NUUR
      print *,"FAKE ADJOINT"
      BS=AS
      RETURN
      END

      SUBROUTINE STOREP(I)
      IMPLICIT NONE
#include "truncation.h"
      INTEGER I
c     print *,"FAKE STOREP"
      RETURN
      END

      SUBROUTINE TMALO3SP(AS,BS)
      IMPLICIT NONE
#include "truncation.h"
      REAL*8 AS(NSH2,3),BS(NSH2,3)
      BS=AS
      RETURN
      END

      SUBROUTINE ALOMA3SP(AS,BS)
      IMPLICIT NONE
#include "truncation.h"
      REAL*8 AS(NSH2,3),BS(NSH2,3)
      INTEGER NUUR
      BS=AS
      RETURN
      END

 
