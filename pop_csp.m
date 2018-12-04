% pop_csp() - Determine optimal common spatial patterns to discriminate
%             between two datasets
%
% Usage:
%   >> [ALLEEG, W]=pop_csp(ALLEEG, datasetlist, chansubset, chansubset2, trainingwindowlength, trainingwindowoffset)
%
%   Inputs:
%       ALLEEG               - array of datasets
%       datasetlist          - list of datasets
%       chansubset           - vector of channel subset for dataset 1          [1:min(Dset1,Dset2)]
%       chansubset2          - vector of channel subset for dataset 2          [chansubset]
%       trainingwindowlength - Length of training window in samples            [all]
%       trainingwindowoffset - Offset(s) of training window(s) in samples      [1]
%
%       Implementation of 
%       Herbert Ramoser, Johanner Muller-Gerking, and Gert Pfurtscheller
%       Optimal Spatial Filtering of Single Trial EEG During Imagined Hand Movement
%       IEEE Transactions on Rehabilitation Engineering
%       Volume 8, Number 4, pp. 441-446, December 2000
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)
%          and Paul Sajda (ps629@columbia,edu 2004)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson, Lucas Parra and Paul Sajda
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [ALLEEG, W]=pop_csp(ALLEEG, setlist, chansubset, chansubset2, trainingwindowlength, trainingwindowoffset)

%%% Preferences %%%
zeromean=1;
%%%%%%%%%%%%%%%%%%%

if nargin < 1
    help pop_csp;
    return;
end;   
if isempty(ALLEEG)
    error('pop_csp(): cannot process empty sets of data');
end
if nargin < 2
    % which set to save
    % -----------------
    promptstr    = { 'Enter two datasets to compare (ex: 1 3):' ...
            'Enter channel subset ([] = all):' ...
            'Enter channel subset for dataset 2 ([] = same as dataset 1):' ...
            'Training Window Length (samples [] = all):' ...
            'Training Window Offset (samples, 1 = epoch start):' };
    inistr       = { '1 2' '' '' '' '1' };
    result       = inputdlg2( promptstr, 'Common Spatial Patterns -- pop_csp()', 1,  inistr, 'pop_csp');
    if length(result) == 0 return; end;
    setlist   	 = eval( [ '[' result{1} ']' ] );
    chansubset   = eval( [ '[' result{2} ']' ] );
    if isempty( chansubset ), chansubset = 1:min(ALLEEG(setlist(1)).nbchan,ALLEEG(setlist(2)).nbchan); end;
    
    chansubset2  = eval( [ '[' result{3} ']' ] );
    if isempty(chansubset2), chansubset2=chansubset; end
    
    trainingwindowlength    = str2num(result{4});
    if isempty(trainingwindowlength)
        trainingwindowlength = ALLEEG(setlist(1)).pnts;
    end;
    trainingwindowoffset    = str2num(result{5});
    if isempty(trainingwindowoffset)
        trainingwindowoffset = 1;
    end;
end;

try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;

if max(trainingwindowlength+trainingwindowoffset-1)>ALLEEG(setlist(1)).pnts,
    error('pop_csp(): training window exceeds length of dataset 1');
end
if max(trainingwindowlength+trainingwindowoffset-1)>ALLEEG(setlist(2)).pnts,
    error('pop_csp(): training window exceeds length of dataset 2');
end
if (length(chansubset)~=length(chansubset2)),
    error('Number of channels from each dataset must be equal.');
end


%%%%%%%%%%%%%%%%%%%

chansubset = 1:min(ALLEEG(setlist(1)).nbchan,ALLEEG(setlist(2)).nbchan);
chansubset2=chansubset;

%%% CSP %%%

% Ramoser equation (1)

t=trainingwindowoffset;
duration=trainingwindowlength;

for i=1:2,

    for trial=1:ALLEEG(i).trials,

        % fprintf('\rProcessing Spatial Covariance - Dataset %d - Tau %d - Trial %d/%d',i,t,trial,ALLEEG(i).trials);
        E=squeeze(ALLEEG(setlist(i)).data(:,t:t+duration-1,trial)); % channels x samples
        [chans,T]=size(E);

        if zeromean, E=E-repmat(mean(E,2),[1 T]); end

        tmpC = (E*E');
        C{i}(:,:,trial) = tmpC./trace(tmpC);

    end % trial
    Cmean{i}=mean(C{i},3);

end % i

% Ramoser equation (2)
Ccomposite=Cmean{1}+Cmean{2};

% Sort eigenvalues in descending order
[Ucomposite,Lambdacomposite] = eig(Ccomposite);
[Lambdacomposite,ind] = sort(diag(Lambdacomposite),'descend');
Ucomposite = Ucomposite(:,ind);

% Ramoser equation (3) - Whitening transform
P=sqrt(inv(diag(Lambdacomposite)))*Ucomposite';

% Ramoser equation (4)
S{1}=P*Cmean{1}*P';
S{2}=P*Cmean{2}*P';

% Ramoser equation (5)
[B,D] = eig(S{1},S{2}); % Simultanous diagonalization
			% Should be equivalent to [B,D]=eig(S{1});

[D,ind] = sort(diag(D)); B = B(:,ind);

W=(B'*P); % Projection matrix

for i=1:length(ind), W(i,:)=W(i,:)./norm(W(i,:)); end

A=pinv(W); % Common spatial patterns

for i=1:2,
    ALLEEG(setlist(i)).icaweights=zeros(ALLEEG(setlist(i)).nbchan);
    
    % In case a subset of channels are used, assign unused electrodes in scalp projection to NaN
    ALLEEG(setlist(i)).icawinv=nan.*ones(ALLEEG(setlist(i)).nbchan);
    ALLEEG(setlist(i)).icasphere=eye(ALLEEG(setlist(i)).nbchan);
end

ALLEEG(setlist(1)).icaweights(chansubset,chansubset)=W;
ALLEEG(setlist(2)).icaweights(chansubset2,chansubset2)=W;
ALLEEG(setlist(1)).icawinv(chansubset,chansubset)=A; 
ALLEEG(setlist(2)).icawinv(chansubset2,chansubset2)=A; 

eeg_options; 
for i=1:2,
    if option_computeica
        ALLEEG(setlist(i)).icaact    = (ALLEEG(setlist(i)).icaweights*ALLEEG(setlist(i)).icasphere)*reshape(ALLEEG(setlist(i)).data, ALLEEG(setlist(i)).nbchan, ALLEEG(setlist(i)).trials*ALLEEG(setlist(i)).pnts);
        ALLEEG(setlist(i)).icaact    = reshape( ALLEEG(setlist(i)).icaact, size(ALLEEG(setlist(i)).icaact,1), ALLEEG(setlist(i)).pnts, ALLEEG(setlist(i)).trials);
    end;
end

for i=1:2, [ALLEEG]=eeg_store(ALLEEG,ALLEEG(setlist(i)),setlist(i)); end

fprintf('Done.\n');
