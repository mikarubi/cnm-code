function [W0, B0, HistoryThreshold, HistoryMeanDelt] = ...
    simann_constraint_model(W0, B0, W, B, D, M, Nl, Nr, objS, objL, objR, IdxN, IdxM, opts)

if ~exist('opts', 'var')
    thr0      = 1;
    thr_decy  = 1-(1e-3);
    thr1      = 0.005;
    time_totl = 1e+9;
    time_decy = 1e+4;
    rng_seed  = rand;
    binary_cn = 1;
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
B0(~W0) = 0;

HistoryMeanDelt  = nan(1,ceil(time_totl/time_decy));
HistoryThreshold = nan(1,ceil(time_totl/time_decy));

n = length(W);
W0(1:n+1:end) = 0;
B0(1:n+1:end) = 0;
W( 1:n+1:end) = 0;
B( 1:n+1:end) = 0;
D( 1:n+1:end) = 0;

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

if ~isequal(sort(W(:)),sort(W0(:)))
    error('Inequality of weights.')
end
if ~exist('IdxN', 'var')
    if objS || objR
        IdxN = true(1, 2*n);
    else
        IdxN = false(1, 2*n);
    end
end
if ~exist('IdxM', 'var')
    switch objL
        case 2;
            IdxM = true(m);
        case 1;
            IdxM = eye(m);
        otherwise;
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
