%% this example shows how to use staircase toolbox
clear all
% initialize it
initializeStaircase
% exp routine
ntrail = 0;
while (1)

    r = showStimuliAndGetResp(mysc.xnext);
    mysc.update(r);
end