function [sigBands, starText, yOffset] = add_significance_to_plot(response, prestim_feature, site, var_method, plots_show, bands_of_interest, reftep)
% Add significance to plots

% 1=theta, 2=alpha, 3=beta, 4=gamma

% STATS CURRENTLY ONLY FOR ALPHA + BETA!!!!!!!

if reftep
    if strcmpi(plots_show, 'VARrel')

        if strcmpi(bands_of_interest, 'alphabeta') || strcmpi(bands_of_interest, 'alpha')

            if strcmpi(response, 'MEP')

                if strcmpi(prestim_feature, 'Power')

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {'**'};

                            sigBands.lowrand = [2];
                            starText.lowrand = {'**'};

                            sigBands.highrand = [2];
                            starText.highrand = {'*'};

                            yOffset = 8;

                        elseif strcmpi(var_method, 'SD')

                            sigBands.lowhigh = [2 3];
                            starText.lowhigh = {'*','**'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [2 3];
                            starText.highrand = {'*','*'};

                            yOffset = 8;

                        elseif strcmpi(var_method, 'IQR')
                            sigBands.lowhigh = [3];
                            starText.lowhigh = {'*'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [3];
                            starText.highrand = {'*'};

                            yOffset = 8;
                        elseif strcmpi(var_method, 'CQV')
                            sigBands.lowhigh = [3];
                            starText.lowhigh = {'**'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [3];
                            starText.highrand = {'*'};

                            yOffset = 8;
                        end
                    elseif strcmp(site, 'Aalto')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {'***'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end

                    elseif strcmp(site, 'combined')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2 3];
                            starText.lowhigh = {'***', '**'};

                            sigBands.lowrand = [2];
                            starText.lowrand = {'*'};

                            sigBands.highrand = [2];
                            starText.highrand = {'**'};

                            yOffset = 8;
                        end
                    end

                elseif strcmpi(prestim_feature, 'Connectivity')

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {'**'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [2];
                            starText.highrand = {'***'};

                            yOffset = 8;

                        elseif strcmpi(var_method, 'SD')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {'**'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [2];
                            starText.highrand = {'***'};

                            yOffset = 8;

                        elseif strcmpi(var_method, 'IQR')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        elseif strcmpi(var_method, 'CQV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    elseif strcmp(site, 'Aalto')
                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    elseif strcmp(site, 'combined')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {'***'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [2];
                            starText.highrand = {'*'};

                            yOffset = 8;
                        end
                    end

                elseif strcmpi(prestim_feature, 'Phase')
                    if strcmpi(site, 'Tuebingen')
                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        elseif strcmpi(var_method, 'SD')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        elseif strcmpi(var_method, 'IQR')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        elseif strcmpi(var_method, 'CQV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    elseif strcmp(site, 'Aalto')
                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    elseif strcmp(site, 'combined')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    end

                elseif strcmpi(prestim_feature, 'PowerConnectivity')

                    % Low Low - High Low
                    % Low Low - Low High
                    % Low Low - High High
                    % High Low - Low High
                    % High Low - High High
                    % Low High - High High

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [3];
                            starText.asterix = {'**'};

                            yOffset = 8;

                        elseif strcmpi(var_method, 'SD')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 8;
                        end

                    elseif strcmp(site, 'Aalto')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 8;
                        end

                    elseif strcmp(site, 'combined')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [1, 3, 6];
                            starText.asterix = {'*', '**', '*'};

                            yOffset = 8;
                        end
                    end

                elseif strcmpi(prestim_feature, 'PowerPhase')

                    % Trough Low - High
                    % Rising Low - High
                    % Peak Low - High
                    % Falling Low - High

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [1 2 3 4];
                            starText.asterix = {'**', '**', '+', '*'};

                            yOffset = 8;

                        elseif strcmpi(var_method, 'SD')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 8;
                        end
                    elseif strcmp(site, 'Aalto')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [1];
                            starText.asterix = {'**'};

                            yOffset = 8;
                        end

                    elseif strcmp(site, 'combined')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [1 2 3 4];
                            starText.asterix = {'***', '**', '*', '**'};

                            yOffset = 8;
                        end
                    end
                end
            end
        end
    end

    if strcmpi(plots_show, 'VARrel')

        if strcmpi(bands_of_interest, 'thetaalphabetagamma')

            if strcmpi(response, 'MEP')

                if strcmpi(prestim_feature, 'Power')

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {"**"};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    end

                elseif strcmpi(prestim_feature, 'Connectivity')

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2];
                            starText.lowhigh = {'***'};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [2];
                            starText.highrand = {'***'};

                            yOffset = 8;
                        end
                    end

                elseif strcmpi(prestim_feature, 'Phase')
                    if strcmpi(site, 'Tuebingen')
                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 8;
                        end
                    end
                elseif strcmpi(prestim_feature, 'PowerConnectivity')

                    % Low Low - High Low
                    % Low Low - Low High
                    % Low Low - High High
                    % High Low - Low High
                    % High Low - High High
                    % Low High - High High

                    if strcmpi(site, 'Tuebingen')
                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 8;
                        end
                    end
                elseif strcmpi(prestim_feature, 'PowerPhase')

                    % Trough Low - High
                    % Rising Low - High
                    % Peak Low - High
                    % Falling Low - High

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 8;
                        end
                    end
                end
            end
        end
    end

    if strcmpi(plots_show, 'MEAN_AMP')

        if strcmpi(bands_of_interest, 'alphabeta') || strcmpi(bands_of_interest, 'alpha')

            if strcmpi(response, 'MEP')

                if strcmpi(prestim_feature, 'Power')

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [2 3];
                            starText.lowhigh = {'***', '**'};

                            sigBands.lowrand = [2];
                            starText.lowrand = {'*'};

                            sigBands.highrand = [2];
                            starText.highrand = {'*'};

                            yOffset = 0.1;
                        end
                    end

                elseif strcmpi(prestim_feature, 'Connectivity')

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.lowhigh = [];
                            starText.lowhigh = {};

                            sigBands.lowrand = [];
                            starText.lowrand = {};

                            sigBands.highrand = [];
                            starText.highrand = {};

                            yOffset = 0.1;
                        end
                    end

                elseif strcmpi(prestim_feature, 'Phase')
                    if strcmpi(site, 'Tuebingen')
                        if strcmpi(var_method, 'CV')
                            sigBands = [1 3; 1 4; 2 3; 3 5];
                            starText = {'*', '*', '*', '*'};

                            yOffset = 0.031;
                        end
                    end
                elseif strcmpi(prestim_feature, 'PowerConnectivity')

                    % Low Low - High Low
                    % Low Low - Low High
                    % Low Low - High High
                    % High Low - Low High
                    % High Low - High High
                    % Low High - High High

                    if strcmpi(site, 'Tuebingen')
                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 0.1;
                        end
                    end
                elseif strcmpi(prestim_feature, 'PowerPhase')

                    % Trough Low - High
                    % Rising Low - High
                    % Peak Low - High
                    % Falling Low - High

                    if strcmpi(site, 'Tuebingen')

                        if strcmpi(var_method, 'CV')
                            sigBands.pair = [];
                            starText.asterix = {};

                            yOffset = 0.1;
                        end
                    end
                end
            end
        end
    end
end

if ~reftep
    if strcmpi(plots_show, 'VARrel')

        if strcmpi(bands_of_interest, 'alphabeta') || strcmpi(bands_of_interest, 'alpha')

            if strcmpi(response, 'MEP')

                if strcmpi(prestim_feature, 'Power')

                    sigBands.lowhigh = [2 3];
                    starText.lowhigh = {'***', '**'};

                    sigBands.lowrand = [2];
                    starText.lowrand = {'**'};

                    sigBands.highrand = [2];
                    starText.highrand = {'*'};

                    yOffset = 8;

                elseif strcmpi(prestim_feature, 'Connectivity')

                    sigBands.lowhigh = [];
                    starText.lowhigh = {};

                    sigBands.lowrand = [];
                    starText.lowrand = {};

                    sigBands.highrand = [];
                    starText.highrand = {};

                    yOffset = 8;

                elseif strcmpi(prestim_feature, 'Phase')

                    sigBands.lowhigh = [];
                    starText.lowhigh = {};

                    sigBands.lowrand = [];
                    starText.lowrand = {};

                    sigBands.highrand = [];
                    starText.highrand = {};

                    yOffset = 8;
                end
            end
        end
    end
end

end


