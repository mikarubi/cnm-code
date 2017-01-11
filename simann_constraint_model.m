function [W0, B0, HistoryThreshold, HistoryMeanDelt] = ...
    simann_constraint_model(W, B, D, M, Nl, Nr, objS, objL, objR, IdxN, IdxM, opts)
%SIMANN_CONSTRAINT_MODEL     Unbiased sampling of networks with hard constraints
%
%   [W0, B0, HistoryThreshold, HistoryMeanDelt] = simann_constraint_model(  ...
%                           W, B, D, M, Nl, Nr, objS, objL, objR, IdxN, IdxM, opts)
%
%   This function randomizes the input network via minimization with
%   simulated weighted) of weighted (+/- binary) node-strength, wiring-cost
%   and module-weight constraints).
%
%   Inputs (for a network with n nodes, m modules and c constraints):
%
%       samp,   Number of networks to sample.
%
%       W,      (length n) square directed and weighted connectivity
%               matrix. All weights must be nonnegative real numbers.
%
%       B,      (length n) square directed and weighted bandwidth matrix.
%               Bandwidth represents an estimate for the cross-section of
%               the white-matter pathway in computations of wiring costs.yes,
%               In the absence of such estimates, the weights matrix W is
%               typically used as a substitute.
%
%       D,      (length n) square directed and weighted distance matrix.
%               This matrix represents an estimate for the distance between
%               all pairs of nodes.
%
%       M,      (length n) module affiliation vector. This vector is often
%               obtained as the output of a community detection algorithm.
%               The vector must contain nonnegative integers, with zeros
%               specifying nodes which are not part of any community. This
%               input may be left empty if there are no module constraints.
%
%       Nl,     (length n) logical vector. This vector denotes regions in
%               the left hemisphere with 1 (and with 0 otherwise).
%
%       Nr,     (length n) logical vector. This vector denotes regions in
%               the right hemisphere with 1 (and with 0 otherwise).
%
%       objS,   flag for constraining node strengths (0 or 1).
%
%       objL,   flag for constraining intra-module weights (0 or 1).
%
%       objR,   flag for constraining wiring-costs (0 or 1).
%
%       IdxN,   (length n) logical vector of node indices to constrain.
%               Alternatively, specify 1 to constrain all nodes or 0 to
%               constrain no nodes. Empty or absent input results in
%               default behavour (no constraints).
%
%       IdxM,	(length m) logical matrix of module indices to constrain.
%               Alternatively, specify 2 to constrain all inter-module and
%               intra-module weights, 1 to constrain all intra-module
%               weights, or 0 for no constraints. Empty or no input results
%               in default behavour (no constraints).
%
%       opts,   custom simulated annealing options 
%               (this argument may be omitted, see code for details).
%
%
%   Outputs:
%
%       W0,     randomized weights matrix.
%
%       B0,     randomized bandwidth matrix.
%
%
%   Notes:
%
%   	The mex file needs to be compiled as: mex loop_sa.c rand_mt.c.
%
%       Nl and Nr, region affiliations for the left and right hemispheres,
%       are specified when one wishes to preserve hemispheric symmetry. One
%       must specify the same number of regions on the left and right, and
%       in the same order. For example, the regions in the matrix may be
%       ordered as follows: L1, L2, ..., Lk, R1, R2, ..., Rk. Regions which
%       are not specified to be part of the left or right hemispheres are
%       not subject to hemispheric symmetry constraints.
%
%
%   Example:
%               % get community structure of a weighted network W
%               M = community_louvain(W, 2);
%
%               % specify node and module constraints
%               n = length(W);                          % number of nodes
%               m = max(M);                             % number of modules
%               Nl = false(n, 1); Nl(1:n/2) = 1;      	% nodes on the left hemisphere
%               Nr = false(n, 1); Nr(1:n/2) = 1;        % nodes on the right hemisphere
%               objS = 1;                               % constrain strengths
%               objL = 1;                               % constrain intra-module weights
%               objR = 1;                               % constrain wiring costs
%               IdxN = true(n, 1);                      % constrain all nodes
%               IdxM = eye(m);                          % constrain all intra-module weights
%
%               % sample networks with the above constraints
%               [W0, B0, HistoryThreshold, HistoryMeanDelt] = ...
%                   simann_constraint_model(  ...
%                           W, B, D, M, Nl, Nr, objS, objL, objR, IdxN, IdxM)
%
%
%   References: Schreiber (1998) Phys Rev Lett 80 2105
%               Rubinov (2016) Nat Commun 7:13812
%
%
%   2016, Mika Rubinov, Janelia HHMI

%   Modification History
%   Dec 2016: Original.


if ~exist('opts', 'var')
    thr0      = 1;                	% initial cutoff error
    thr_decy  = 1-(1e-3);          	% temperature decay rate
    thr1      = 0.005;             	% target cutoff error
    time_totl = 1e+9;              	% total number of iterations
    time_decy = 1e+4;              	% number of iterations at each temperature
    rng_seed  = rand;              	% initial seed
    binary_cn = 1;                 	% binary constraints (on or off)
else
    thr0        = opts.thr0;
    thr_decy    = opts.thr_decy;
    thr1        = opts.thr1;
    time_totl   = opts.time_totl;
    time_decy   = opts.time_decy;
    rng_seed    = opts.rng_seed;
    binary_cn   = opts.binary_cn;
end

B(  ~W) = 0;

HistoryMeanDelt  = nan(1,ceil(time_totl/time_decy));
HistoryThreshold = nan(1,ceil(time_totl/time_decy));

n = length(W);
W( 1:n+1:end) = 0;
B( 1:n+1:end) = 0;
D( 1:n+1:end) = 0;
W0 = 1 * W + 0;
B0 = 1 * B + 0;


M = M(:);
m = max(M);

[I, J]  = find(~eye(n));

Ns     = (Nr | Nl);                        % symmetric nodes
Na     = ~Ns;
Ms     = zeros(size(Ns));                  % symmetry map
Ms(Nl) = find(Nr);
Ms(Nr) = find(Nl);
INa    = find(Na);                         % indices of asymmetric nodes
n_asym = numel(INa);
n_edge = numel(I);
n_obj  = nnz(isfinite(nonzeros([objS objL objR])));

if ~exist('IdxN', 'var')
    if objS || objR
        IdxN = true(1, 2*n);
    else
        IdxN = false(1, 2*n);
    end
end
if ~exist('IdxM', 'var')
    switch objL
        case 2
            IdxM = true(m);
        case 1
            IdxM = eye(m);
        otherwise
            IdxM = false(m);
    end
end

IdxN      = int32(IdxN);
IdxM      = int32(IdxM);
M_c       = int32(M-1);                     % convert to C indexing
binary_cn = int32(logical(binary_cn));
objS      = int32(objS);
objL      = int32(objL);
objR      = int32(objR);
n_obj     = int32(n_obj);
Ns        = int32(Ns);
Na        = int32(Na);
Ms_c      = int32(Ms-1);                    % convert to C indexing
INa_c     = int32(INa-1);                   % convert to C indexing
I_c       = int32(I-1);                     % convert to C indexing
J_c       = int32(J-1);                     % convert to C indexing
n         = int32(n);
n_edge    = int32(n_edge);
m         = int32(m);
n_asym    = int32(n_asym);
time_totl = int32(time_totl);
time_decy = int32(time_decy);
rng_seed  = int32(rng_seed);

temp      = double(1*thr0);                 % for SA only
temp_decy = double(1*thr_decy);             % for SA only

loop_sa;
