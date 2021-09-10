import("stdfaust.lib");


// -------------------------------------------------------------------
// This file demonstrates how to use the Faust FDS library to couple
// multiple Finite Difference Schemes at arbitrary points with rigid
// connections. The overall logic is to define an enviroment like
// coupledSchemes, to be found at the end of the file, which contains
// all individual FDS models constructed by the FDS library as well as
// the coupling information.
//
// `system1D` is then the method that performs the update and coupling
// equations, making use of `forceUpdate` to calculate and route the 
// forces after performing the individual update steps of the FDS.
//
// The example in this file demonstrates the coupling of three stiff 
// strings over a bridge, but the code was kept as general as possible
// to allow arbitrary couplings.
// -------------------------------------------------------------------

declare name "Coupled Finite Difference Schemes in Faust";
declare version "0.1";
declare author "David Suedholt";

// the samplerate here is hardcoded to calculate nPoints from L at compile time
k = 1.0 / 44100;

// basic stencil parameters 
nNeighbors = 2; // R
nTimesteps = 1; // T


// this is fd.model1D without the recursion, so that we can add the forces before going to the next time step
// A coupled system consists of one schemeUpdate per FDS stacked in parallel
schemeUpdate(points,R,T,scheme) = fd.route1D(points,R,T,scheme) : fd.buildScheme1D(points,R,T);


//---------------------------`forceUpdate1D`---------------------------------------
// Given a number of FDS schemes whose individual update equations have already
// been calculated, calculate the force at each coupling point and add it back 
// to the affected grid points according to the order of the coupling
// (positive force to the system coupled 'above' the other, negative force to
// the one coupled 'below' the other)
//
// #### Usage
//
// ```
// si.bus(totalPoints) : forceUpdate1D(coupledSchemes) : si.bus(totalPoints);
// ```
//
// Where:
//
// * `coupledSchemes`: An environment containing the information about the schemes
// and their coupling as defined further below
//------------------------------------------------------------------------------
forceUpdate(coupledSchemes) =
    inRouting : interp : forcecalc : outRouting
with {
    M = coupledSchemes.nSchemes;
    nPoints = coupledSchemes.nPoints;
    startingPoint = coupledSchemes.startingPoint;
    totalPoints = startingPoint(M); // sum of all points because of the definition of startingPoint
    h = coupledSchemes.h;
    posR = coupledSchemes.posR;
    posS = coupledSchemes.posS;

    beta = coupledSchemes.beta;

    nCouplings = coupledSchemes.nCouplings;

    r = coupledSchemes.r;
    s = coupledSchemes.s;

    // alpha values for the interpolation and spreading operators
    alphaR(i) = posR(i) - floor(posR(i));
    alphaS(i) = posS(i) - floor(posS(i));

    // the force calculation depends on two grid points from each scheme, so four points per force
    // this routing attaches these points to the bottom of the scheme so that they can be interpolated and used further
    inRouting = route(totalPoints, totalPoints+4*nCouplings,
        // simply pass the current scheme through
        par(i, totalPoints, (i+1, i+1)),

        // attach all relevant points from all couplings
        par(i, nCouplings, 
            // the two points from the 'above' part of the coupling, i.e the r-th scheme
            (startingPoint(r(i)) + floor(posR(i)), totalPoints+i*4+1), (startingPoint(r(i)) + floor(posR(i))+1, totalPoints+i*4+2),
            // the two points from the 'below' part, i.e. the s-th scheme
            (startingPoint(s(i)) + floor(posS(i)), totalPoints+i*4+3), (startingPoint(s(i)) + floor(posS(i))+1, totalPoints+i*4+4)   
        )
    );

    // linearly interpolate each pair of grid points, so that we are down to one value per scheme
    interp = 
        // still just passing the scheme through
        si.bus(totalPoints),
        
        par(i, nCouplings, 
            // (1 - alpha) * u_l + (alpha) * u_{l+1}
            +(*(1 - alphaR(i)), *(alphaR(i))), +(*(1 - alphaS(i)), *(alphaS(i)))
        );

    // now for the actual force calculation
    forcecalc =
        // attach the beta coefficients for all couplings
        (si.bus(totalPoints+2*nCouplings), par(i, nCouplings, beta(i))) :
        
        // still just passing the scheme through
        si.bus(totalPoints),

        // route each coefficient to its interpolated grid point, perform the multiplication
        // and add each two grid points belonging to one coupling to obtain the forces
        (ro.interleave(2*nCouplings, 2) : par(i, 2*nCouplings, *) : par(i, nCouplings, +) <:

        // perform the spreading operation: each force is distributed back to the four grid points it affects
        ro.interleave(nCouplings, 4) : par(i, nCouplings,
            // the force is added to the two points in the r-th scheme
            *((1 - alphaR(i))/h(r(i))),*(alphaR(i)/h(r(i))), 
            // and subtracted from the ones in the s-th scheme
            *((alphaS(i) - 1)/h(s(i))),*((alphaS(i) * -1)/h(s(i))))
        );


    // finally the calculated forces are routed back to the grid points they affect
    outRouting = route(totalPoints+4*nCouplings, totalPoints,
        // the original scheme, to which the force is added by the routing
        par(i, totalPoints, (i+1, i+1)),

        // route the force back to the grid
        par(i, nCouplings,
            (totalPoints+i*4+1, startingPoint(r(i)) + floor(posR(i))), (totalPoints+i*4+2, startingPoint(r(i)) + floor(posR(i))+1),
            (totalPoints+i*4+3, startingPoint(s(i)) + floor(posS(i))), (totalPoints+i*4+4, startingPoint(s(i)) + floor(posS(i))+1)
        )
    );
};



//---------------------------`system1D`---------------------------------------
// Given a number of FDS schemes stacked in parallel and their coupling points,
// calculate the individual updates and add the coupling forces at each time step.
// Takes as input an external force signal (e.g. excitation) for each point in the 
// combined scheme.
//
// #### Usage
//
// ```
// si.bus(totalPoints) : system1D(coupledSchemes) : si.bus(totalPoints);
// ```
//
// Where:
//
// * `coupledSchemes`: An environment containing the information about the schemes
// and their coupling as defined further below
//------------------------------------------------------------------------------
system1D(coupledSchemes) =
    (routing : schemes : forceUpdate(coupledSchemes) : norm) ~ si.bus(totalPoints)
with {
    M = coupledSchemes.nSchemes;
    nPoints = coupledSchemes.nPoints;
    startingPoint = coupledSchemes.startingPoint;
    totalPoints = startingPoint(M);
    schemes = coupledSchemes.schemes;
    a = coupledSchemes.a;

    // after adding all forces, divide by the supplied factor
    norm = par(i, M, par(j, nPoints(i), /(a(i))));

    // the inputs are coming as points0, points1, ..., forces0, forces1 ...
    // this routing rearranges them to points0, forces0, points1, forces1 ...
    routing = route(totalPoints*2, totalPoints*2,
        // routing of grid points
        par(i, M, par(j, nPoints(i), (startingPoint(i)+j+1, startingPoint(i)*2+j+1))),
        // routing of forces
        par(i, M, par(j, nPoints(i), (startingPoint(i)+totalPoints+j+1, startingPoint(i)*2+nPoints(i)+j+1)))
    );
};


// Coefficients for the FDS of a stiff string with the simply supported boundary condition
// In the simply supported condition, we define the coefficients at the boundary points such that
// u_{-1} = -u_3 and u_{N+1} = -u_{N-3}
simplySupportedScheme(params) = coeffsLeft, midCoeffsDelay, par(i, nPoints-2, midCoeffs, midCoeffsDelay), coeffsRight, midCoeffsDelay
with {
    nPoints = params.nPoints;
    EI = params.E * params.I;
    h = params.h;
    rhoA = params.rho * params.Area;
    s0 = params.sigma0;
    s1 = params.sigma1;
    T = params.T;

    coeffsLeft =    0,
                    0,
                    (-2*T*h^2 - 6*EI - 2*rhoA*s0*h^4/k - 4*rhoA*s1*h^2/k + 2*rhoA*h^4/k^2 + EI) / h^4,
                    (T*h^2 + 4*EI + 2*rhoA*s1*h^2/k) / h^4,
                    (-EI) / h^4;

    midCoeffs = (-EI) / h^4,
                (T*h^2 + 4*EI + 2*rhoA*s1*h^2/k) / h^4,
                (-2*T*h^2 - 6*EI - 2*rhoA*s0*h^4/k - 4*rhoA*s1*h^2/k + 2*rhoA*h^4/k^2) / h^4,
                (T*h^2 + 4*EI + 2*rhoA*s1*h^2/k) / h^4,
                (-EI) / h^4;

    coeffsRight =   (-EI) / h^4,
                    (T*h^2 + 4*EI + 2*rhoA*s1*h^2/k) / h^4,
                    (-2*T*h^2 - 6*EI - 2*rhoA*s0*h^4/k - 4*rhoA*s1*h^2/k + 2*rhoA*h^4/k^2 + EI) / h^4,
                    0,
                    0;
    
    midCoeffsDelay = 0,
                     (-2*rhoA*s1/h^2/k), 
                     (2*rhoA*s0/k + 4*rhoA*s1/h^2/k - rhoA/k^2),
                     (-2*rhoA*s1/h^2/k),
                     0;
};

// Coefficients for the FDS of a stiff string with the clamped supported boundary condition
// In the clamped condition, no special treatment for the boundary points is needed, the routing
// just supplies the virtual points as zero
clampedScheme(params) = par(i, nPoints, coeffs, coeffsDelay)
with {
    nPoints = params.nPoints;
    EI = params.E * params.I;
    h = params.h;
    rhoA = params.rho * params.Area;
    s0 = params.sigma0;
    s1 = params.sigma1;
    T = params.T;

    coeffs = (-EI) / h^4,
                (T*h^2 + 4*EI + 2*rhoA*s1*h^2/k) / h^4,
                (-2*T*h^2 - 6*EI - 2*rhoA*s0*h^4/k - 4*rhoA*s1*h^2/k + 2*rhoA*h^4/k^2) / h^4,
                (T*h^2 + 4*EI + 2*rhoA*s1*h^2/k) / h^4,
                (-EI) / h^4;
    
    coeffsDelay = 0,
                     (-2*rhoA*s1/h^2/k), 
                     (2*rhoA*s0/k + 4*rhoA*s1/h^2/k - rhoA/k^2),
                     (-2*rhoA*s1/h^2/k),
                     0;
};

// This is a convenience function to calculate all required parameters based on the "tunable" parameter set further down
calcAllParams(params) = environment {
    // keep all parameters from the initial structure
    rho = params.rho;
    E = params.E;
    L = params.L;
    f0 = params.f0;
    radius = params.radius;
    sigma0 = params.sigma0;
    sigma1 = params.sigma1;
    I = params.I;
    Area = params.Area;

    T = (2*f0*L)^2 * rho * Area; // Tension based on the desired pitch

    // Calculate the minimum grid spacing given by the stability condition
    hmin = sqrt(k/2*(T*k/rho/Area + 4*sigma1 + sqrt((T*k/rho/Area + 4*sigma1)^2 + 16*E*I/rho/Area)));
    
    // nPoints is the number of points that are actually updated each time step.
    // N = floor(L / hmin) assumes grid points u_0, u_1, ..., u_{N-1}, u_N 
    // where u_0 = u_N = 0 are fixed boundary points, so we are left with nPoints = N-1 "updatable" points
    nPoints = floor(L / hmin) - 1; 
    // Calculate the actual grid spacing based on the number of points
    h = L / (nPoints+1);
};

string1 = environment {
    rho = 1200; // Material density
    E = 2e9;    // Young's modulus
    L = 0.65;   // Length of the string
    radius = 4.6e-4; 
    f0 = 146.83;     // Fundamental frequency, used for tension calculation   
    sigma0 = 1.38;   // Frequency-independent damping coefficient
    sigma1 = 1.3e-4; // Frequency-dependent damping coefficient

    I = ma.PI * radius^4 / 4; // Moment of Inertia for the strings
    Area = ma.PI * radius^2; // Cross-sectional area for the strings
};

string2 = environment {
    rho = 1200;       
    E = 2e9;       
    L = 0.65;         
    radius = 4.6e-4;      
    f0 = 169;        
    sigma0 = 1.38;      
    sigma1 = 1.3e-4;    
    I = ma.PI * radius^4 / 4; 
    Area = ma.PI * radius^2;
};

string3 = environment {
    rho = 1200;       
    E = 2e9;       
    L = 0.65;         
    radius = 4.6e-4;      
    f0 = 246.94;      
    sigma0 = 1.38;      
    sigma1 = 1.3e-4;    
    I = ma.PI * radius^4 / 4; 
    Area = ma.PI * radius^2;
};

bridge = environment {
    rho = 1500;
    E = 3e9;
    I = 1.136e-10; // Moment of inertia for the bridge
    Area = 2e-5; // cross-sectional area for the bridge
    L = 0.16;         
    radius = 4.6e-4; 
    sigma0 = 1.343;
    sigma1 = 2.59e-3;
    f0 = 0; // results in Tension being set to 0
};

params(0) = calcAllParams(string1);
params(1) = calcAllParams(string2);
params(2) = calcAllParams(string3);
params(3) = calcAllParams(bridge);


// Calculation of the beta coefficients needed for force calculation
calcBeta(r, s, posr, poss, ar, as) = (-1)/denom/ar, 1/denom/as
with {
    hr = params(r).h;
    hs = params(s).h;

    alphaR = posr - floor(posr);
    alphaS = poss - floor(poss);

    jnormsqr = ((1 - alphaR)^2 + alphaR^2) / hr^2;
    jnormsqs = ((1 - alphaS)^2 + alphaS^2) / hs^2;

    denom = hr * jnormsqr / ar + hs * jnormsqs / as;
};

// Define the coupling of the r-th scheme at position x_r above the s-th scheme at position x_s
// in this case, the three strings (0 - 2) are coupled at x = 10 cm 
// to the bridge (3), where they are fixed at 4cm, 8cm and 12 cm
couplings = environment {
    nCouplings = 3;
    
    // indices of the coupled systems (coupling r above s)
    r(0) = 0;
    s(0) = 3;

    // positions are in meters
    xr(0) = 0.1; 
    xs(0) = 0.04;

    r(1) = 1;
    s(1) = 3;

    xr(1) = 0.1;
    xs(1) = 0.08;

    r(2) = 2;
    s(2) = 3;

    xr(2) = 0.1;
    xs(2) = 0.12;
};

// The "master" environment containing all the FDS schemes and the coupling information
coupledSchemes = environment {
    nSchemes = 4; // total number of schemes

    // Number of points in each scheme
    nPoints(i) = params(i).nPoints;
    
    // stack all schemes in parallel using functions from the fd library
    schemes = par(i, nSchemes-1, schemeUpdate(params(i).nPoints, nNeighbors, nTimesteps, simplySupportedScheme(params(i)))),
                schemeUpdate(params(3).nPoints, nNeighbors, nTimesteps, clampedScheme(params(3)));

    // since the points of all schemes are stacked together we need to 
    // determine where one schemes stops and the next begins
    startingPoint(0) = 0;
    startingPoint(i) = startingPoint(i-1) + nPoints(i-1);

    // the normalization factors; division by a is the last step in every update
    a(i) = params(i).rho * params(i).Area / k^2;
    
    // grid spacing is needed for force calculation and spreading operators
    h(i) = params(i).h;

    // copy over the coupling definitions
    nCouplings = couplings.nCouplings;
    r = couplings.r;
    s = couplings.s;
    xr = couplings.xr;
    xs = couplings.xs;

    // calculate the continuous position in points based on the position in meters
    posR(i) = xr(i) / h(r(i));
    posS(i) = xs(i) / h(s(i));

    // the tuple beta_r, beta_s of coefficients for the force calculation
    beta(i) = calcBeta(r(i), s(i), posR(i), posS(i), a(r(i)), a(s(i)));
};


//----------------------------Interface Elements-----------------------------//

play(i) = button("Play%i");
inPoint(i) = hslider("Input Point%i",floor(params(i).nPoints/2),0,params(i).nPoints-1,1);
outPoint(i) = hslider("Output Point%i",floor(params(i).nPoints/2),0,params(i).nPoints-1,0.01):si.smoo;
outSlider = hslider("Select output to listen", 0, 0, coupledSchemes.nSchemes-1, 1);
maxforce = hslider("Maximal Excitation Force", 1, 0.1, 10, 0.1);

//------------------------------Excitation Force-----------------------------//

excdur_sec = 0.005;
excdur_samp = int(excdur_sec * 44100);
excforce(n) = ba.if(n <= excdur_samp, maxforce/2 * (1 - cos(ma.PI * n/excdur_samp)), 0);

// ba.countup will continue to put out the upper limit once it is reached, until triggered again
// so have the limit be one higher and check in excforce if we've exceeded it
forceModel(i) = (ba.countup(excdur_samp+1, 1-(play(i))):excforce)/params(i).h;

//----------------------------------Process---------------------------------//

input = par(i, coupledSchemes.nSchemes, (forceModel(i)/params(i).h)<:fd.linInterp1D(params(i).nPoints, inPoint(i)));
output = par(i, coupledSchemes.nSchemes, fd.linInterp1DOut(params(i).nPoints, outPoint(i)));
outSelect = ba.selectn(coupledSchemes.nSchemes, outSlider);

process =  input : system1D(coupledSchemes) : output : outSelect <: _,_;
