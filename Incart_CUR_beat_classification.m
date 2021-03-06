% This script applies 1 nearest neighbor classification to the CUR beat
% annotation results on the Incart Database to
% classify those beats not selected by CUR.

% Results are generated for the original PhysioNet annotations, as well as
% for the AAMI and AAMI2 labels. 

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create cell array of the different record names
patient_ID_list = cell(75,1);

for id = 1:75
    if id > 9
        patient_ID_list{id} = ['I' num2str(id) 'm'];
    else
        patient_ID_list{id} = ['I0' num2str(id) 'm'];
    end
end

% The CUR tolerances tested in beat selection (recall that these values
% are divided by 10 prior to being input into the incremental QR code)
CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];

% Classify unlabled beats for each patient record and each CUR tolerance
full_classification = cell(1,length(CUR_stopping_tol));
full_true_class = cell(1,length(CUR_stopping_tol));

for tol = 1:length(CUR_stopping_tol)
    
    full_classification{tol} = {};
    full_true_class{tol} = {};
    
    for i = 1:length(patient_ID_list)
        
        load(['Incart_database/' patient_ID_list{i} '_filtered_data_matrix'])

    data_matrix = data_matrix_beat_normalization(info.data_matrix{2}); % looking at Lead II
    
    beat_annotations = info.annotations;
    
    load([patient_ID_list{i} '_CUR_annotation_tracking'])
        
        
        % Perform 1-NN classification, using annotations from the
        % CUR-selected beats to the label the rest of the data set
        k = 1;
        distance = 'euclidean';
        to_match = data_matrix;
        library = data_matrix(annotations.CUR_q{1,tol},:);
        [k_nn{i,tol},k_smallest_distances] = find_k_nearest_neighbors(to_match,library,k,distance);
        
        % Assign label according to nearest CUR-selected beat
        classification{i,tol} = beat_annotations(annotations.CUR_q{1,tol}(k_nn{i,tol}));
        
        full_classification{tol} = [full_classification{tol} classification{i,tol}];
        full_true_class{tol} = [full_true_class{tol} beat_annotations];
        
    end
    
    
    load('Incart_CUR_correct_annotation_ident_results.mat')
    
    % Form a confusion matrix to evaluate classification results 
    C{tol} = confusionmat(full_true_class{tol},full_classification{tol},'order',results.total_annotations);
    
    % Add column and row totals to confusion matrix
    C_summary{tol} = C{tol};
    C_summary{tol}(:,end+1) = sum(C{tol},2);
    C_summary{tol}(end+1,:) = sum(C_summary{tol},1);
    
    % Compute sensitivity and positive predictive value for each label
    S_i(tol,:) = diag(C{tol})./C_summary{tol}(1:end-1,end)*100; % sensitivity
    P_i(tol,:) = (diag(C{tol})')./C_summary{tol}(end,1:end-1)*100; % positive predictive value
    
    % Compute global sensitivity, positive predictive value, and accuracy 
    S(tol) = nansum(diag(C{tol})./C_summary{tol}(1:end-1,end))/length(results.total_annotations)*100; % global sensitivity
    P(tol) = nansum((diag(C{tol})')./C_summary{tol}(end,1:end-1))/length(results.total_annotations)*100; % global positive predictive value
    A(tol) = sum(diag(C{tol}'))/C_summary{tol}(end,end)*100; % global accuracy
    
    %% Translate PhysioNet classification results to AAMI and AAMI2 labling schemes
    % AAMI Classes:
    % % N = N, L, R, e, j
    % % S = A, a, J, S
    % % V = V, E
    % % F = F
    % % Q = /, f, Q
    
    % AAMI2 Classes:
    % % N = N, L, R, e, j
    % % S = A, a, J, S
    % % V = V, E, F
    
    % % Q = /, f, Q (Discarded for MIT database...)
    
    AAMI_annotations = {'N' 'S' 'V' 'F' 'Q'};
    AAMI2_annotations = {'N' 'S' 'V_hat' 'Q'};
    
    % Map true labels and predicted labels to AAMI and AAMI2 annotations
    AAMI_full_true_class{tol} = cell(size(full_true_class{tol}));
    AAMI2_full_true_class{tol} = cell(size(full_true_class{tol}));
    
    AAMI_full_classification{tol} = cell(size(full_classification{tol}));
    AAMI2_full_classification{tol} = cell(size(full_classification{tol}));
    for r = 1:length(full_true_class{tol})
        if strcmp(full_true_class{tol}{r},'N') || strcmp(full_true_class{tol}{r},'L') || strcmp(full_true_class{tol}{r},'R') || strcmp(full_true_class{tol}{r},'e') || strcmp(full_true_class{tol}{r},'j') || strcmp(full_true_class{tol}{r},'B') || strcmp(full_true_class{tol}{r},'n')
            AAMI_full_true_class{tol}{1,r} = AAMI_annotations{1};
            AAMI2_full_true_class{tol}{1,r} = AAMI2_annotations{1};
        elseif strcmp(full_true_class{tol}{r},'A') || strcmp(full_true_class{tol}{r},'a') || strcmp(full_true_class{tol}{r},'J') || strcmp(full_true_class{tol}{r},'S')
            AAMI_full_true_class{tol}{1,r} = AAMI_annotations{2};
            AAMI2_full_true_class{tol}{1,r} = AAMI2_annotations{2};
        elseif strcmp(full_true_class{tol}{r},'V') || strcmp(full_true_class{tol}{r},'E')
            AAMI_full_true_class{tol}{1,r} = AAMI_annotations{3};
            AAMI2_full_true_class{tol}{1,r} = AAMI2_annotations{3};
        elseif strcmp(full_true_class{tol}{r},'F')
            AAMI_full_true_class{tol}{1,r} = AAMI_annotations{4};
            AAMI2_full_true_class{tol}{1,r} = AAMI2_annotations{3};
        elseif strcmp(full_true_class{tol}{r},'/') || strcmp(full_true_class{tol}{r},'f') || strcmp(full_true_class{tol}{r},'Q')
            AAMI_full_true_class{tol}{1,r} = AAMI_annotations{5};
            AAMI2_full_true_class{tol}{1,r} = AAMI2_annotations{4};
        elseif strcmp(full_true_class{tol}{r},'!')
            AAMI_full_true_class{tol}{1,r} = AAMI_annotations{6};
            AAMI2_full_true_class{tol}{1,r} = AAMI2_annotations{5};
        end
        
        if strcmp(full_classification{tol}{r},'N') || strcmp(full_classification{tol}{r},'L') || strcmp(full_classification{tol}{r},'R') || strcmp(full_classification{tol}{r},'e') || strcmp(full_classification{tol}{r},'j') || strcmp(full_classification{tol}{r},'B') || strcmp(full_classification{tol}{r},'n')
            AAMI_full_classification{tol}{1,r} = AAMI_annotations{1};
            AAMI2_full_classification{tol}{1,r} = AAMI2_annotations{1};
        elseif strcmp(full_classification{tol}{r},'A') || strcmp(full_classification{tol}{r},'a') || strcmp(full_classification{tol}{r},'J') || strcmp(full_classification{tol}{r},'S')
            AAMI_full_classification{tol}{1,r} = AAMI_annotations{2};
            AAMI2_full_classification{tol}{1,r} = AAMI2_annotations{2};
        elseif strcmp(full_classification{tol}{r},'V') || strcmp(full_classification{tol}{r},'E')
            AAMI_full_classification{tol}{1,r} = AAMI_annotations{3};
            AAMI2_full_classification{tol}{1,r} = AAMI2_annotations{3};
        elseif strcmp(full_classification{tol}{r},'F')
            AAMI_full_classification{tol}{1,r} = AAMI_annotations{4};
            AAMI2_full_classification{tol}{1,r} = AAMI2_annotations{3};
        elseif strcmp(full_classification{tol}{r},'/') || strcmp(full_classification{tol}{r},'f') || strcmp(full_classification{tol}{r},'Q')
            AAMI_full_classification{tol}{1,r} = AAMI_annotations{5};
            AAMI2_full_classification{tol}{1,r} = AAMI2_annotations{4};
        elseif strcmp(full_classification{tol}{r},'!')
            AAMI_full_classification{tol}{1,r} = AAMI_annotations{6};
            AAMI2_full_classification{tol}{1,r} = AAMI2_annotations{5};
        end
    end

    % Construct AAMI results confusion matrix
    AAMI_C{tol} = confusionmat(AAMI_full_true_class{tol},AAMI_full_classification{tol},'order',AAMI_annotations);
    % Add column and row totals to confusion matrix
    AAMI_C_summary{tol} = AAMI_C{tol};
    AAMI_C_summary{tol}(:,end+1) = sum(AAMI_C{tol},2);
    AAMI_C_summary{tol}(end+1,:) = sum(AAMI_C_summary{tol},1);
    
    % Compute sensitivity and positive predictive value for each label
    AAMI_S_i(tol,:) = diag(AAMI_C{tol})./AAMI_C_summary{tol}(1:end-1,end)*100; % sensitivity
    AAMI_P_i(tol,:) = (diag(AAMI_C{tol})')./AAMI_C_summary{tol}(end,1:end-1)*100; % positive predictive value
    
    % Compute global sensitivity, positive predictive value, and accuracy
    AAMI_S(tol) = nansum(diag(AAMI_C{tol})./AAMI_C_summary{tol}(1:end-1,end))/length(AAMI_annotations)*100; % global sensitivity
    AAMI_P(tol) = nansum((diag(AAMI_C{tol})')./AAMI_C_summary{tol}(end,1:end-1))/length(AAMI_annotations)*100; % global postive predictive value
    AAMI_A(tol) = sum(diag(AAMI_C{tol}'))/AAMI_C_summary{tol}(end,end)*100; % global accuracy
    
    
    % Construct AAMI2 results confusion matrix
    AAMI2_C{tol} = confusionmat(AAMI2_full_true_class{tol},AAMI2_full_classification{tol},'order',AAMI2_annotations);
    % Add column and row totals to confusion matrix
    AAMI2_C_summary{tol} = AAMI2_C{tol};
    AAMI2_C_summary{tol}(:,end+1) = sum(AAMI2_C{tol},2);
    AAMI2_C_summary{tol}(end+1,:) = sum(AAMI2_C_summary{tol},1);
    
    % Compute sensitivity and positive predictive value for each label
    AAMI2_S_i(tol,:) = diag(AAMI2_C{tol})./AAMI2_C_summary{tol}(1:end-1,end)*100; % sensitivity
    AAMI2_P_i(tol,:) = (diag(AAMI2_C{tol})')./AAMI2_C_summary{tol}(end,1:end-1)*100; % positive predictive value
    
    % Compute global sensitivity, positive predictive value, and accuracy 
    AAMI2_S(tol) = nansum(diag(AAMI2_C{tol})./AAMI2_C_summary{tol}(1:end-1,end))/length(AAMI2_annotations)*100; % global sensitivity
    AAMI2_P(tol) = nansum((diag(AAMI2_C{tol})')./AAMI2_C_summary{tol}(end,1:end-1))/length(AAMI2_annotations)*100; % global postive predictive value
    AAMI2_A(tol) = sum(diag(AAMI2_C{tol}'))/AAMI2_C_summary{tol}(end,end)*100; % global accuracy
    
    % Generate positive predictive value, sensitivity, and F1 score on V
    % and V_hat labels
    
    P_V(tol) = AAMI_P_i(tol,3);
%     P_V(tol) = AAMI_C_summary{tol}(3,3)/AAMI_C_summary{tol}(end,3)*100;
    
    S_V(tol) = AAMI_S_i(tol,3);
%     S_V(tol) = AAMI_C_summary{tol}(3,3)/AAMI_C_summary{tol}(3,end)*100;
    
    P_V_hat(tol) = AAMI2_P_i(tol,3);
%     P_V_hat(tol) = AAMI2_C_summary{tol}(3,3)/AAMI2_C_summary{tol}(end,3)*100;
    
    S_V_hat(tol) = AAMI2_S_i(tol,3);
%     S_V_hat(tol) = AAMI2_C_summary{tol}(3,3)/AAMI2_C_summary{tol}(3,end)*100;
    
    AAMI_F_V(tol) = 2*AAMI_C_summary{tol}(3,3)/(AAMI_C_summary{tol}(3,end)+AAMI_C_summary{tol}(end,3))*100;
    
    AAMI2_F_V_hat(tol) = 2*AAMI2_C_summary{tol}(3,3)/(AAMI2_C_summary{tol}(3,end)+AAMI2_C_summary{tol}(end,3))*100;
end

% Output results
full_summary = [S' P' A'];
AAMI_summary = [AAMI_S' AAMI_P' AAMI_A'];
AAMI2_summary = [AAMI2_S' AAMI2_P' AAMI2_A'];

tolerance_labels = {'0.05' '0.01' '5 \times 10^{-3}' '1 \times 10^{-3}' '5 \times 10^{-4}' '1 \times 10^{-4}' '5 \times 10^{-5}' '1 \times 10^{-5}'}';

andper = repmat('&',length(CUR_stopping_tol),1);
endline = repmat('\\',length(CUR_stopping_tol),1);
% PhysioNet classification results
table_full_results = table(tolerance_labels,andper,S',andper,P',andper,A',andper,S_i(:,1),andper,P_i(:,1),andper,S_i(:,2),andper,P_i(:,2),andper,S_i(:,3),andper,P_i(:,3),andper,S_i(:,4),andper,P_i(:,4),andper,S_i(:,5),andper,P_i(:,5),andper,S_i(:,6),andper,P_i(:,6),andper,S_i(:,7),andper,P_i(:,7),andper,S_i(:,8),andper,P_i(:,8),andper,S_i(:,9),andper,P_i(:,9),andper,S_i(:,10),andper,P_i(:,10),endline);
table_full_results.Properties.VariableNames = {'tol' 'andper1' 'Se' 'andper2' 'P_plus' 'andper' 'A' 'andper3' ['Se_' results.total_annotations{1}] 'andper4' ['P_' results.total_annotations{1}] 'andper5' ['Se_' results.total_annotations{2}] 'andper6' ['P_' results.total_annotations{2}] 'andper7' ['Se_' results.total_annotations{3}] 'andper8' ['P_' results.total_annotations{3}] 'andper9' ['Se_' results.total_annotations{4}] 'andper10' ['P_' results.total_annotations{4}] 'andper11' ['Se_' results.total_annotations{5}] 'andper12' ['P_' results.total_annotations{5}] 'andper13' ['Se_' results.total_annotations{6}] 'andper14' ['P_' results.total_annotations{6}] 'andper15' ['Se_' results.total_annotations{7}] 'andper16' ['P_' results.total_annotations{7}] 'andper17' ['Se_' results.total_annotations{8}] 'andper18' ['P_' results.total_annotations{8}] 'andper19' ['Se_' results.total_annotations{9}] 'andper20' ['P_' results.total_annotations{9}] 'andper21' ['Se_' results.total_annotations{10}] 'andper22' ['P_' results.total_annotations{10}] 'endline'};
writetable(table_full_results,'Incart_full_CUR_classification_tol_tests','Delimiter',' ','WriteRowNames',true)

% AAMI classification results
table_AAMI_results = table(tolerance_labels,andper,AAMI_S',andper,AAMI_P',andper,AAMI_A',andper,AAMI_S_i(:,1),andper,AAMI_P_i(:,1),andper,AAMI_S_i(:,2),andper,AAMI_P_i(:,2),andper,AAMI_S_i(:,3),andper,AAMI_P_i(:,3),andper,AAMI_S_i(:,4),andper,AAMI_P_i(:,4),andper,AAMI_S_i(:,5),andper,AAMI_P_i(:,5),endline);
table_AAMI_results.Properties.VariableNames = {'tol' 'andper' 'Se' 'andper1' 'P_plus' 'andper0' 'A' 'andper2' ['Se_' AAMI_annotations{1}] 'andper3' ['P_' AAMI_annotations{1}] 'andper4' ['Se_' AAMI_annotations{2}] 'andper5' ['P_' AAMI_annotations{2}] 'andper6' ['Se_' AAMI_annotations{3}] 'andper7' ['P_' AAMI_annotations{3}] 'andper8' ['Se_' AAMI_annotations{4}] 'andper9' ['P_' AAMI_annotations{4}] 'andper10' ['Se_' AAMI_annotations{5}] 'andper11' ['P_' AAMI_annotations{5}] 'endline'};
writetable(table_AAMI_results,'Incart_AAMI_CUR_classification_tol_tests','Delimiter',' ','WriteRowNames',true)

% AAMI2 classification results
table_AAMI2_results = table(tolerance_labels,andper,AAMI2_S',andper,AAMI2_P',andper,AAMI2_A',andper,AAMI2_S_i(:,1),andper,AAMI2_P_i(:,1),andper,AAMI2_S_i(:,2),andper,AAMI2_P_i(:,2),andper,AAMI2_S_i(:,3),andper,AAMI2_P_i(:,3),andper,AAMI2_S_i(:,4),andper,AAMI2_P_i(:,4),endline);
table_AAMI2_results.Properties.VariableNames = {'tol' 'andper' 'Se' 'andper1' 'P_plus' 'andper0' 'A' 'andper2' ['Se_' AAMI2_annotations{1}] 'andper3' ['P_' AAMI2_annotations{1}] 'andper4' ['Se_' AAMI2_annotations{2}] 'andper5' ['P_' AAMI2_annotations{2}] 'andper6' ['Se_' AAMI2_annotations{3}] 'andper7' ['P_' AAMI2_annotations{3}] 'andper8' ['Se_' AAMI2_annotations{4}] 'andper9' ['P_' AAMI2_annotations{4}] 'endline'};
writetable(table_AAMI2_results,'Incart_AAMI2_CUR_classification_tol_tests','Delimiter',' ','WriteRowNames',true)

% V and V_hat results for parameter selection
table_AAMI_V_hat_and_V_results = table(tolerance_labels,andper,AAMI2_S_i(:,3),andper,AAMI2_P_i(:,3),andper,AAMI2_F_V_hat',andper,AAMI_S_i(:,3),andper,AAMI_P_i(:,3),andper,AAMI_F_V',endline);
table_AAMI_V_hat_and_V_results.Properties.VariableNames = {'tol' 'andper1' 'Se_Vhat' 'andper2' 'P_plus_Vhat' 'andper3' 'F_1_Vhat' 'andper4' 'Se_V' 'andper5' 'P_plus_V' 'andper6' 'F_1_V' 'endline'};
writetable(table_AAMI_V_hat_and_V_results,'Incart_AAMI_V_hat_and_V_CUR_classification_tol_tests','Delimiter',' ','WriteRowNames',true)

% Save confusion matrices
save('Incart_confusion_matrices','C_summary','AAMI_C_summary','AAMI2_C_summary')


