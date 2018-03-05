%This code sets up parameter structures for a better set of dynamic experiments. This code sets up a parameter array as required by matlabs parfor loop. This code calls the integrating function and parses the data.


%This is the structure that gets passed to the actual solver. This is a "blank copy" that gets copied and filled in for each simulation. These are the values that are constant for ALL simulations, and a few that are easily generated from web data.
params = struct(... 
    ... Web Parameters:
    'S',[]...
    ,'C',[]...
    ... BS Parameters:
    ,'K',[]...
    ... FR Parameters:
    ,'halfSat',[]...
    ,'phi',[]...
    ,'h',[]...
    ,'eijBas',[] ...
    ,'eijCon',[] ...
    ,'yFree',[] ...
    ,'yPara',[] ...
    ,'axFree',[]...
    ,'axPara',[]...
    ... Web Properties:
    ... Link Level properties::
    ,'res',[]...
    ,'con',[]...
    ,'eij',[]...
    ,'wij',[] ...
    ,'yij',[]...
    ... Species level Properties::
    ,'basal',zeros(40,1)...
    ,'r',zeros(40,1)...
    ,'B0',zeros(40,1)...
    ,'x',zeros(40,1)...
    ,'M',zeros(40,1)...
    ,'modelCode',[0 0]...
    ,'para',zeros(40,1)...
    ... Simulation Parameter:
    ,'Tf',[]...
    ,'Trelax',[]...
    ... Solver Parameters:
    ,'extctThresh',[]...
    ,'AbsTol',[]...
    ,'RelTol',[]...
    ,'odeSolver',[]... .
    ,'options',[]...
    );
%Web Parameters
S = 40;
C = 0.15;
nBasal = 6;

params.S = S;
params.C = C;

%BS Parameters
K = 5;

params.K = K;

%FR Parameters
halfSat = 0.5;
phi = 0.15;
h = 1.2;
eijBas = 0.45;
eijCon = 0.85;
yFree = 8;
yPara = yFree;
axFree = 0.314;
axPara = axFree;

params.halfSat = halfSat;
params.phi = phi;
params.h = h;
params.eijBas = eijBas;
params.eijCon = eijCon;
params.yFree = yFree;
params.yPara = yPara;
params.axFree = axFree;
params.axPara = axPara;

%Simulation Parameters: 
TFinal = 10000;
TRelax =  2000;

params.Tf = TFinal;
params.Trelax = TRelax;

%Solver Parameters;
extctThresh = 1e-10;
AbsTol = 1e-16;
RelTol = 1e-6;
options = odeset()

params.extctThresh = extctThresh;
paramsAbsTol = 1e-16;
params.RelTol = 1e-6;
params.options = options;

%Factors
nWeb = 100;
useBioModel = [false true];
fPar = linspace(1,17,9)/34;

%This cell array holds holds the data that describe the levels of all factors tat are being tested in this simulation. 
simParams = cell(nSims,1);

[simParams{:}] = deal(struct(...
                 'web',0 ...
                ,'fPar',0 ...
                ,'useBioModel',0 ...
                ,'LL',[] ...
                ,'B0',zeros(40,1) ...
                ,'gr',zeros(40,1) ...
                ,'parOrder', zeros(40,1) ...
                ,'patl', zeros(40,1)
                ));

                
