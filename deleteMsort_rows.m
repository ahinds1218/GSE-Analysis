%Deleting the rows
numParticles = max(M2(:,3))+1
n = numParticles;
total_frames = max(M2(:,4))+1
k = total_frames;

M2((((k).*(n))+1):end,:) = [];

M = M2;


