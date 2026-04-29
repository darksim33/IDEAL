function [FitResults,FitQuality,Params] = ideal_fitter(Data_raw, Mask_raw, Data_raw_masked, Params)


%% Perform IDEAL fitting
tStart = tic;
for slice = Params.slice

    % Select single slice for image resizing
    Data = squeeze(Data_raw(:,:,slice,:));
    Data_masked = squeeze(Data_raw_masked(:,:,slice,:));
    Mask = squeeze(Mask_raw(:,:,slice,:));

    % prepare Output structs
    FitQuality = cell(size(Data_raw,3),1,1);
    ROIstat = cell(size(Data_raw,3),1,1);

    % Check if sclice contains data after masking. If not skip slice
    if isempty(Data_masked(~isnan(Data_masked)))
        FitQuality{slice} = {};
        ROIstat{slice} = {};
        fprintf('Slice %s does not contain any data after masking and will be skipped.\n',...
            num2str(slice));
        continue
    end

    % Loop Resampling steps
    for res_step = 1:size(Params.Dims_steps,1)
        fprintf('Downsampling step no: %s of slice %s\n', num2str(res_step), num2str(slice));

        if res_step == 1 && Params.mean
        % Average in first step
        Mask_res = ones(1,1);
        Data_res = zeros(1, 1, numel(Params.b_values));

        for i = 1:numel(Params.b_values)
                Data_res(:,:,i) = nanmean(Data_masked(:,:,i),'all');
        end

        % Resample Data if not in finale step
        else if res_step < size(Params.Dims_steps, 1)
                [Mask_res, Data_res] = ...
                    resample_data(Mask,Data, ...
                    [Params.Dims_steps(res_step,1), Params.Dims_steps(res_step,2)],...
                    Params.b_values);
            else
                Data_res = Data_masked;
                Mask_res = Mask;
            end
        end
        % Prepare fitting parameters
        [fit, FitResults, gof, output, op] = setup_fitting(Params,res_step,size(Mask_res));

        % Start Voxelvise-Fitting
        for x = 1:size(Mask_res,1)
            for y = 1:size(Mask_res,2)

                % if matrix size is not 1x1, update parameters
                if res_step ~= 1
                    op = update_fitting(op, fit_res, Params, x, y);
                end
                if Mask_res(x,y)
                    if ~isnan(Data_res(x,y,1))
                        % only run if voxel is part of mask and data
                        % is present else skip
                        switch Params.Model
                            case {"biexp","Biexp"}
                                [FitResults{x,y}, gof{x,y}, ...
                                    output{x,y}] = ...
                                    BiexpFit(Params.b_values, ...
                                    squeeze(Data_res(x,y,:))', op);
                            case {"biexp_T1corr","Biexp_T1corr"}
                                [FitResults{x,y}, gof{x,y}, ...
                                    output{x,y}] = ...
                                    BiexpFit_T1corr(Params.b_values, ...
                                    squeeze(Data_res(x,y,:))', op, Params.TM, Params.TR, Params.TE);
                            case {"triexp","Triexp"}
                                [FitResults{x,y}, gof{x,y}, ...
                                    output{x,y}] = ...
                                    TriexpFit(Params.b_values, ...
                                    squeeze(Data_res(x,y,:))', op);
                        end
                        fit = extract_fit(fit,FitResults, ...
                                        Params.Model,x,y);

                    end
                end
            end
        end
        % interpolate parameters
        if res_step < size(Params.Dims_steps,1)
            fit_res = interpolate_fit(fit,Params.Dims_steps, ...
                                        res_step,Params.Model);
        end
    end
    fprintf("Fitting Completed!\nStarting Plotting...\n");
    Params.time = toc(tStart);
