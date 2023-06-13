%   ---------------------------------------------------------------
%   Function Name:  main

function main()

% redeclarations of globals
global data 
global IG
IG=[1 2 3 4 5 6]';

global metofcal          % metofcal=1
global N_pop             % N_pop=50%input('no pop:');
global N_HL              % N_HL=2%input('Hidden Layer:');
global pmm               % pmm=.1%input('Mutation prob:'); 
global N_Generation      % N_Generation=150%input('Generation:') 
global N_Prediction      % N_Prediction=5%input('Prediction:'); 
global first_seed        % first_seed=-400%input('seed:');
global CR                % Crossover constant
global F                 % Scaling factor

ReadData();

% Call main function
[y,bestfitness,final]=main12(data,N_pop,N_HL,pmm,N_Generation,N_Prediction,metofcal,first_seed,CR,F);

