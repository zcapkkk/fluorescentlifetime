function [output]=groupproject_simulate(amplitudes,lifetimes,acquisitiontime)
    bg=10; %number of background counts per second, keep at 10
    
    %check that each amplitude has a corresponding lifetime
    if length(amplitudes)~=length(lifetimes)
        return
    end
    
    %create empty vector to store decay data
    puredecay=zeros(381,1);
    
    %normalise amplitudes, just in case they didn't initially sum to 1
    amplitudes=amplitudes/sum(amplitudes);
    
    %generate a multiexponential decay starting at 1 at t=0
    %using the supplied amplitudes and lifetimes
    for i=1:381
        t=(1/19)*(i-1); %each bin is (1/19) ns, starting at t=0
        for j=1:length(amplitudes)
            puredecay(i,1)=puredecay(i,1)+amplitudes(j)*exp(-t/lifetimes(j));
        end
    end
    
    %we do our measurements at 2500 counts per second
    %calculate how many fluorescence counts per second this corresponds to
    %i.e. subtract background from total counts
    fluorate=2500-bg;
    
    %calculate total number of fluorescence photons counted in measurement
    totalfluorescence=fluorate*acquisitiontime; 
    
    %now scale the multiexponential decay so it contains this many counts
    puredecay=totalfluorescence*puredecay/sum(puredecay);

    %and add on 'bg' counts per second spread evenly across all bins
    puredecay=puredecay+(bg*acquisitiontime/381);
    
    %finally add poisson noise to each bin
    noisydecay=poissrnd(puredecay);
    
    %and tidy up output with a time axis
    output=zeros(381,2);
    for i=1:381
        output(i,1)=(i-1)*(1/19);
        output(i,2)=noisydecay(i,1);
    end
end