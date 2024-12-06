## Fortran codes for Computation of the land-sea unified complete Bouguer effect on gravity outside geoid
https://www.zcyphygeodesy.com/en/h-nd-126.html
## [Algorithm purpose]
    Using the numerical integral method, from the land-sea terrain model and ellipsoidal height grid file of the land-sea surface, compute the land-sea unified complete Bouguer effect on the gravity (mGal) on the geoid or in near-Earth space.
    The program is suitable for the unified computation of the complete Bouguer effect on gravity, gravity anomaly and gravity disturbance in land, land-sea junction and sea area.
    The land terrain complete Bouguer effect is defined as the effect of the terrain mass above the geoid on the Earth's gravity field. While the seawater complete Bouguer effect is defined as the effect on the Earth's gravity field because of the density of seawater compensated to the density of land terrain.
The land-sea unified complete Bouguer effects can also be calculated by spherical harmonic analysis and synthesis of global land-sea terrain masses.
## [Main program for test entrance]
    landseacmpBougeffintgrl.f90
    The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m)......
    Input files: The land-sea terrain model file, in which the land terrain height is greater than zero while the seafloor water depth is less than zero. The ellipsoidal height grid file of the land-sea surface, which stands for the land-sea surface position employed to calculate the integral distance.
    The record format of the output file reslt.txt: Behind the record of the calculation point file, appends the local terrain effect, spherical shell Bouguer effect and sear-water complete Bouguer effect and landsea complete Bouguer effect.
## (1) Algorithm module for numerical integral of local terrain effects on various field elements
    LTerAllBLH(BLH,dtm,sfh,nlat,nlon,hd,dr,GRS,ter)
    Input parameters:BLH(3) - longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m) of the calculation point.
    Input parameters: dtm(nlat,nlon) - the ground digital elevation model (normal /orthometric height) grid (>0), which is employed to indicate land terrain relief.
    Input parameters: sfh(nlat,nlon) - the ground ellipsoidal height grid, which represents the terrian surface position employed to calculate the integral distance.
    Input parameters: dr, hd(6) - the integral radius (m) and grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid).
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Return parameters: ter(5) - local terrain effects on the height anomaly (m), gravity (anomaly/disturbance, mGal), vertical deflection (ʺ, to south, to west) or (disturbing) gravity gradient (E, radial).
## (2) Algorithm module for numerical integral of sear-water complete Bouguer effects
    OwBougreBLH(BLH,sea,sfh,nlat,nlon,hd,dr,GRS,ter)
    Input parameters: sea(nlat,nlon) - the sea-floor depth model grid, which is employed to indicate the sea-floor relief.
    Return parameters: ter(2) - sear-water complete Bouguer effects on the height anomaly (m) and gravity (anomaly/disturbance, mGal).
## (3) Calculation module for the normal gravity field
    normdjn(GRS,djn); GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m²/s²), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (6) Algorithm library for interpolation point value from numerical grid
    CGrdPntD(lon,lat,dt,row,col,hd)；CGrdPntD2(lon,lat,dt,row,col,hd)
    CShepard(lon,lat,dt,row,col,hd)；Gauss2D(lon,lat,dt,row,col,hd)
## (7) Other auxiliary modules
    PickRecord(str0, kln, rec, nn)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    1.4.1 Format convention for geodetic data file
    7.5.1 Expression of land terrain mass gravitational field - land complete Bouguer effect
    7.5.2 Integral formula of local terrain effect outside the Earth
    7.6.1 Marine terrain gravitational field - seawater complete Bouguer effect
    7.1(4) Low-dgree Legendre function and its first and second derivative algorithms

DOS executable test file and all input and output data.
