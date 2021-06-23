c - Dimension of the ft12_diffint files
    
c       IMPLICIT NONE
       INTEGER NX,NY,NZ,IMH
       INTEGER NNT, NXH, NXP2, NYP

       PARAMETER ( NNT=128 )

C       PARAMETER ( NX=256, NY=128, NZ=129 )
       PARAMETER ( NX=512, NY=256, NZ=257 )

       PARAMETER (NXH=NX/2+1,NXP2=NX+2)
       PARAMETER (NYP = NY/NNT)

