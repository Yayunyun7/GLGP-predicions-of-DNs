function [aloc, bloc, eloc, nloc, p2loc, NonDN, DN, DNloc] = PPG_abcde_detection2(ppg, R, Q, Fs)

tic
ppg = ppg(:);

if length(find(isnan(ppg)))
    fprintf('There are NAN in the signal. Correct by setting them to 0...\n') ;
    ppg(isnan(ppg)) = 0;
end

switch nargin
    case 1
        error('Specify the sampling rate Fs as the second argument.')
end

% parameters
threshold = 0.02;
W1len = ceil(0.11*Fs) ;
W2len = ceil(0.66*Fs) ;

W3len = ceil(0.05*Fs) ;
W4len = ceil(0.22*Fs) ;




% filter signals and get APG(accelerated PPG, second derivatives)
[b, a] = butter(2, [0.5, 8] / (Fs / 2));
filtered = filtfilt(b, a, ppg);
APG = Fs * Fs * diff(filtered, 2) ;
APG = [APG; APG(end); APG(end)] ;






%% find a (the same as Q in our notation)
ddfiltered1 = APG ;
ddfiltered1(find(APG<0)) = 0 ;
u = ddfiltered1.^2 ;

V = filtfilt(ones(1, W2len) ./ W2len, 1, u);
indicator = filtfilt(ones(1, W1len) ./ W1len, 1, u);

% setup dynamic threshold
mu = smooth(u, Fs*20, 'loess') ;

try threshold = filtfilt(ones(1, 5 * Fs) ./ (5 * Fs), 1, u) * threshold;
catch
    threshold = threshold .* mu ; %mean(u);
end
% setup box function (ROI)
t = indicator > V + threshold;

% get all supports
[M, start] = regexp(sprintf('%i', [0 t']), '1+', 'match');

M = cellfun(@length, M);


%

aloc = nan(length(R), 1);
adpM = zeros(size(R)) ;
for i = 1:length(R)
    tmp = start - R(i) ; tmp(find(tmp>0)) = inf ; 
    [~, tmpi] = min(abs(tmp)) ;
    %[~, tmpi] = min(abs(start - R(i))) ;
    adpM(i) = M(tmpi) ;
    
    inds = max(1, R(i) - adpM(i)*2 + 1) ;
    [~, aloc(i)] = max(APG(inds:R(i)));
    aloc(i) = inds - 1 + aloc(i) ;

    if i > 1 & aloc(i) == aloc(i-1)
       aloc(i) = aloc(i) + 1 ; 
    end
end








%% find b

bloc = aloc ;
for jj = 1:length(bloc)
    kk = ceil(0.008*Fs)+1 ;
    while kk < R(jj) - aloc(jj) %ceil(0.136*Fs)
        if abs(APG(aloc(jj)+kk))>abs(APG(aloc(jj)+kk-1)) & ...
                abs(APG(aloc(jj)+kk))>abs(APG(aloc(jj)+kk+1))
            break;
        else
            kk = kk + 1 ;
        end
    end
    bloc(jj) = aloc(jj) + kk ;

    if jj > 1 & bloc(jj) == bloc(jj-1)
       bloc(jj) = bloc(jj) + 1 ; 
    end
end






%% find e
% remove a-b
u2 = APG ;
for jj = 1: length(aloc)
    tmp = max(1, Q(jj) - adpM(jj)): min(length(APG), bloc(jj) + round(adpM(jj)/2)) ;
    u2(tmp) = 0;
end
u2(find(u2<0)) = 0  ;
u2 = u2.^2 ;

eloc = nan(length(R), 1);
for i = 1:length(R)
    
    inds = min(length(ppg), R(i) + 1) : min(length(ppg), R(i) + adpM(i)*2 + 1) ;
    [~, tmp] = max(u2(inds));
    eloc(i) = inds(1) - 1 + tmp(1) ;
    
    if i > 1 & eloc(i) == eloc(i-1)
       eloc(i) = eloc(i) + 1 ; 
    end
end








%% find notch & 2nd peak
nloc = eloc ;
p2loc = eloc ;

dppg = diff(filtered) ;
dppg = [dppg; dppg(end)] ;
NonDN = [] ;
DN = [];
DNloc =[];

for jj = 1:length(nloc)
   
    if jj < length(nloc)
        tmp = dppg(eloc(jj)+1: min(aloc(jj+1), eloc(jj)+round(adpM(jj)/2))) ;    
    else
        tmp = dppg(eloc(jj)+1: min(length(dppg), eloc(jj)+round(adpM(jj)/2))) ;
    end
    
    % if there is a notch around e point (if the 1st derivative has a zero)
    % set -0.01 to avoid numerical error
    if length(find(tmp>-0.01))
        
        kk = 1 ;
        
        % find DN
        while  kk<adpM(jj) && nloc(jj)+kk < length(ppg)
            % if found the local minimum, stop.
            if  ppg(nloc(jj)+kk)<ppg(nloc(jj)+kk-2) && ppg(nloc(jj)+kk)<ppg(nloc(jj)+kk+2)
                break;
            else
                kk = kk + 1 ;
            end
        end
        % if the while loop stopped before the possible window ends%
        % which means we find the DN
        if kk<adpM(jj)& nloc(jj)+kk < length(ppg)
            
            % in this case, assign the DN point****
            nloc(jj) = nloc(jj) + kk ;
            DN = [DN jj];
            DNloc = [DNloc nloc(jj)];
            
            %%%%% and determine the 2nd peak
            hh = 1 ;
            while hh < adpM(jj) & nloc(jj)+hh < length(ppg)
                if ppg(nloc(jj)+hh)>ppg(nloc(jj)+hh-2) & ppg(nloc(jj)+hh)>ppg(nloc(jj)+hh+2)
                    break;
                else
                    hh = hh + 1 ;
                end
            end

                      
            % assign the 2nd peak
            p2loc(jj) = nloc(jj) + hh ;
            
            if jj > 1 && nloc(jj) == nloc(jj-1)
                nloc(jj) = nloc(jj) + 1 ; 
            end
            
            if jj > 1 && p2loc(jj) == p2loc(jj-1)
                p2loc(jj) = nloc(jj) + 1 ; 
            end
            
        else
            
            % if cannot find a locally minimal point
            NonDN = [NonDN jj] ;
        end
        
    else  
        % if there is no notch around the e point, record it.
        NonDN = [NonDN jj] ;
    end
end





disp(['a-b-c-d-e detection completed in ' num2str(toc) ' seconds.'])





end


