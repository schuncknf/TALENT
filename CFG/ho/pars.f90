MODULE pars
    integer, parameter:: dp=8
    integer, parameter:: dp2=8
    INTEGER, parameter:: NMAX=100!max n number
    INTEGER, parameter:: LMAX=0!max l number
    INTEGER, parameter:: NGAUSS=100!number of gauss integral points
    REAL(DP),PARAMETER:: PI4=ATAN(1.0D0)
    REAL(DP)::HBOMG,BOSC,HBARM,MASS,HBAR
    REAL(DP)::FACTOR
END MODULE pars
