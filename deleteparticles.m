function [M_sort_cut] = deleteparticles (M_sort)
evalResponse = input('Type the 2-digit particle IDs to be deleted from analysis in a single row matrix, ascending order, each ID separated by one space (Format: [## ## ## ##]).\n');
M_sort_cut = M_sort;
numParticles = max(M_sort(:,3))+1
deletedParticleIDs = evalResponse
numdeletedParticles = length(deletedParticleIDs)    

for k = 1:numdeletedParticles
   M_sort_cut((M_sort(:,3) == evalResponse(k)-k+1),:) = [];
   delindex = deletedParticleIDs(k)+1;
   delindex
   
    for m = delindex:M_sort_cut(end,3)+k-1;
        M_sort_cut((M_sort(:,3)== m-k),3) = m-k; 
    end;
end;

fprintf('Remaining particles:\n')
numParticles = max(M_sort_cut(:,3))+1
pause();

end