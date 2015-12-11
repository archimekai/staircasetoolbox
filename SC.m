classdef SC <handle
    %SC staircase method toolbox
    %   
    %   by WEN Kai    2015-12-09   bilewhale@163.com
    %   Department of Psychology, Peking University
    %   
    %   highlights:  built in graph support
    %                linear or logit mode support
    %                thrshold computation support
    %   ver 0.5
    
    properties
        xnext           % the signal strength for the next trial  if xnext
        ndown           % 
        nup             % 
        step            % a number if mode == 1     a vector if other modes
        steps           % record needed steps   in most conditions need two
        nstep           % optional; how many steps is given in linear condition
        stepratio       % logit stepratio = nextsti / currentsti in down condition
        sti0            % initial sti
        stdsti          % standard sti
        mode            % 0 : linear, no change during whole process   1 : ratio
        ncurtrial       % No. of current Trial + 1   that means the true exp starts from trial 2
        upordowns       % record upordwons for next trial    1: unchanged   2:up  3:down  -1: for trial 1
        reversals       % record this trial is or not a reversal point
        % nreversal       % the number of reversals
        par             % the data structure containing the configurations of the parameter space   
        
        history = struct('resp',[],'sti',[]);         % previous reactions  -1 for trial no. 1
        
        % user data
        userdata01
        userdata02
        userdata03
        userdata04
        
    end
    
    methods
        %constructor
        function sc = SC(par)
            reset(sc, par);
        end
        
        %% INITIALIZATION
        % set parameters
        function reset(sc, par)
            % ?????????? especially the first trial
            sc.par = par;
            sc.xnext = par.sti0;
            sc.sti0 = par.sti0;
            sc.stdsti = par.stdsti;
            sc.nup = par.nup;
            sc.ndown = par.ndown;
            sc.mode = par.mode;
            if(sc.mode == 0)
                sc.step = par.step;
            elseif(sc.mode == 1)
                sc.stepratio = par.stepratio;
            else
                error('wrong mode!');
            end
            
            sc.ncurtrial = 2;  %
            sc.upordowns(1) = -1;
            sc.reversals(1) = 0;
            sc.history.resp(1) = -1;
            sc.history.sti(1) = -1;
            sc.history.sti(2) = sc.sti0;
        end
        
        %% UPDATING REACTION
        function update(sc,r)
            % record reaction
            % r = 1 currect
            sc.history.resp(sc.ncurtrial) = r;
            
            % compute xnext
            if (r == 0) % wrong  may need weaker
                if (min(countSuc(sc, r), nearestUporDown(sc))) >= sc.nup
                    sc.upordowns(sc.ncurtrial) = 2;
                    sc.xnext = nextSti(sc,'up');
                else
                    sc.upordowns(sc.ncurtrial) = 1;
                    sc.xnext = sc.history.sti(sc.ncurtrial);
                end
            elseif (r == 1)
                if (min(countSuc(sc, r), nearestUporDown(sc))) >= sc.ndown
                    sc.upordowns(sc.ncurtrial) = 3;
                    sc.xnext = nextSti(sc,'down');
                    
                    % may be need to control the range of xnext not less
                    % than standrad sti
                else
                    sc.upordowns(sc.ncurtrial) = 1;
                    sc.xnext = sc.history.sti(sc.ncurtrial);
                end
            else
                warning('FATAL ERROR: r must be 0 or 1!');
            end
            sc.reversals(sc.ncurtrial) = isReversal(sc);
            sc.ncurtrial = sc.ncurtrial + 1;
            sc.history.sti(sc.ncurtrial) = sc.xnext;
        end
        
        %% other functions
        function i = isReversal(sc)  % determine if current trial is a reversal 
            % find before trends

                n = 1;
                while(sc.upordowns(end - n) == 1) 
                    n = n + 1;
                end
                laststate = sc.upordowns(end - n);
                if(laststate == -1)
                    i = 0;
                    return;
                else
                    if(sc.upordowns(sc.ncurtrial) == laststate || sc.upordowns(sc.ncurtrial) == 1)
                        i = 0;
                        return;
                    else
                        i = 1;
                        return;
                    end
                end
        end
        function i = countReversal(sc)
            i = sum(sc.reversals);
        end
        %% some private functions
        function i = countSuc(sc,r)  % ?????????????????
            i = 1;
            while(1)
                if (sc.history.resp(sc.ncurtrial - i) == r)
                    i = i + 1;
                else
                    break;
                end
            end                        
        end
%         function i = nearestRev(sc)   %  find the trials between the nearest reversal(not included)  and current trial(included)
%             
%         end
        function i = nearestUporDown(sc)  % find the nearest upordwon point (which indicates its next trial is weaker or stonger than it)
            i = 1;
            while(1)
                if (sc.upordowns(end - i + 1) == 1)
                    i = i + 1;
                else
                    break;
                end
            end
        end
        function i = nextSti(sc,ori)
            if(sc.mode == 0)  % linear, constant step
                if (strcmp(ori,'up'))
                    i = sc.step + sc.history.sti(sc.ncurtrial);
                elseif (strcmp(ori,'down'))
                    i = -1 * sc.step + sc.history.sti(sc.ncurtrial);
                else
                    error('invalid nextstep argument!');
                end
                if(i < sc.stdsti)
                    i = sc.stdsti;
                    warning('too small next sti, reset to standard sti!');
                end
            elseif (sc.mode == 1) % logit  need to know this is up or down,in varargin{1}
                if (strcmp(ori,'up'))
                    i = 1 / sc.stepratio * (sc.history.sti(sc.ncurtrial) - sc.stdsti) + sc.stdsti;
                elseif (strcmp(ori,'down'))
                    i = sc.stepratio * (sc.history.sti(sc.ncurtrial) - sc.stdsti) + sc.stdsti;
                else
                    error('invalid nextstep argument!');
                    % under construction
                end
            end
        end
        
        % compute threshold based on the reversals needed
        % but how to compute?
        % use the mean of the reverse points
        function thr = getthr(sc, reversals)
            revpoints = sc.history.sti(sc.reversals);
            thr = mean(revpoints(end - reversals + 1: end));    
        end
        %% plot a graph
        function plot(sc)
           figure;
           hold on;
           for i = 2:sc.ncurtrial - 1
               if (sc.history.resp(i) == 1)
                    scatter(i, sc.history.sti(i),'ok');
               elseif (sc.history.resp(i) == 0)
                   scatter(i, sc.history.sti(i),'ro');
               else
                   warning('invalid sc.history.sti!');
               end
               if(sc.reversals(i) == 1)
                   scatter(i,sc.history.sti(i),'+b');
               end
           end
           thr = getthr(sc, 6);
           line([thr,thr],[0 0]);
           xlabel('No. of trials');
           ylabel('compare sti');
        end
    end    
end




