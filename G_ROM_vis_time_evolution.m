textFileName = ['GROMQIvis' num2str(IC) '.txt'];
fileID = fopen(textFileName, 'w');  

%initial condition

    ts  = 1;
% 
    b = 1.0/dt * MassROM * velPrevPrev;
    % % build matrix
    NLmat = 0*MassROM;
    for k=1:N
       NLmat = NLmat + velPrevPrev(k)*TriLinROM2(:,:,k);
    end
    
    A = 1.0/dt * MassROM + nu1*StiffROM +nu2 * StiffROMsin +nu3 * StiffROMcos ...
        +nu4 * StiffROMsin2 +nu5 * StiffROMcos2 + NLmat;     
    % 
    A(1,:)=0;
    A(1,1)=1;
    b(1)=1;
    %     
    % % solve the linear system
    velSoln = A \ b;
 
 
   
    vt = 1/dt*(velSoln - velPrevPrev);
    
    lift = -20*( vt'*vlmass' + nu1*velSoln'*vlstiff'+nu2*velSoln'*vlstiffsin' +nu3*velSoln'*vlstiffcos' ...
           + nu4*velSoln'*vlstiffsin2'+nu5*velSoln'*vlstiffcos2'+ velSoln' * NLlift * velSoln );
       
    drag = -20*( vt'*vdmass' + nu1*velSoln'*vdstiff'+ nu2*velSoln'*vdstiffsin' +nu3*velSoln'*vdstiffcos' ...
           + nu4*velSoln'*vdstiffsin2'+ nu5*velSoln'*vdstiffcos2'+ velSoln' * NLdrag * velSoln );


    energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
     
    fprintf(fileID,'%f %f %f %f\n',ts*dt, lift, drag,   energy);

    
    velPrev = velSoln;

    %BDF2 with linear scheme
   
    for ts=2:numTimeSteps
        %RHS
         b = 2.0/dt * MassROM * velPrev -0.5/dt*MassROM*velPrevPrev;
         % build matrix
         NLmat = 0*MassROM;
         for k=1:N
             NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*TriLinROM2(:,:,k);
         end
         
         A = 1.5/dt * MassROM + nu1*StiffROM +nu2 * StiffROMsin +nu3 * StiffROMcos ...
             +nu4 * StiffROMsin2 +nu5 * StiffROMcos2 + NLmat;

         A(1,:)=0;
         A(1,1)=1;
         b(1)=1;
    
         % solve the linear system
         velSoln = A \ b;

         
  
         vt = 1/dt*(velSoln - velPrev);
         
         lift = -20*( vt'*vlmass' + nu1*velSoln'*vlstiff'+nu2*velSoln'*vlstiffsin' ...
                +nu3*velSoln'*vlstiffcos' +nu4*velSoln'*vlstiffsin2'+nu5*velSoln'*vlstiffcos2' ...
                + velSoln' * NLlift * velSoln );
            
         drag = -20*( vt'*vdmass' + nu1*velSoln'*vdstiff'+ nu2*velSoln'*vdstiffsin' ...
                 +nu3*velSoln'*vdstiffcos' + nu4*velSoln'*vdstiffsin2'+ nu5*velSoln'*vdstiffcos2' ...
                 + velSoln' * NLdrag * velSoln );

         energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
     
       
          fprintf(fileID,'%f %f %f %f\n',ts*dt, lift, drag,   energy);
          velPrevPrev=velPrev;
          velPrev = velSoln;
    end