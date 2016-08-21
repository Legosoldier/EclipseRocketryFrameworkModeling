function VAF = setup_aero(instruct)

% Generates important vehicle parameters needed for simulink calculations
% Variables are passed in a struct called VAF, which should be assigned to
% a VAF variable in workspace, which then the VAF_lib should be able to
% access and perform necessary aerodynamic calculations with. 

VAF = {};

switch instruct

    case 'create'
        %% Airframe Parts and Properties
        
        % Nosecone properties
        VAF.NC.LEAD_EDGE = 0;
        VAF.NC.TRAIL_EDGE = 24;
        VAF.NC.TRAIL_EDGE_DIAM = 6.15;
        NC_CREATE('Tangent Ogive', VAF.NC.TRAIL_EDGE_DIAM, VAF.NC.TRAIL_EDGE - VAF.NC.LEAD_EDGE);
        
        % Airframe Properties
        VAF.AF.BODY_DIAM = 6.15;
        VAF.AF.LENGTH = 84;
        VAF.AF.LEAD_EDGE = 24; 
        VAF.AF.TRAIL_EDGE = VAF.AF.LEAD_EDGE + VAF.AF.LENGTH;
        VAF.AF.LEAD_EDGE_AREA = VAF.AF.BODY_DIAM ^2 * 0.25 * pi;
        VAF.AF.TRAIL_EDGE_AREA = VAF.AF.LEAD_EDGE_AREA;
        VAF.AF.X_CP = BODY_CP_CALC(VAF.AF.LENGTH, VAF.AF.TRAIL_EDGE_AREA, VAF.AF.TRAIL_EDGE_AREA*VAF.AF.LENGTH, VAF.AF.LEAD_EDGE_AREA);
        VAF.AF.AREA_PLANAR = VAF.AF.BODY_DIAM * VAF.AF.LENGTH;
        VAF.AF.BODY_RADIUS_AVG = 0.5*(VAF.NC.AREA_PLAN + VAF.AF.AREA_PLANAR) / (VAF.NC.TRAIL_EDGE + VAF.AF.LENGTH);
        VAF.AF.BODY_AREA_WET = VAF.NC.AREA_WET + pi*VAF.AF.BODY_DIAM*VAF.AF.LENGTH;
%         VAF.AF. = (VAF.AF_1.len * VAF.AF_1.A_le)  ;

        % Fin Input Properties
        VAF.FIN.BODY_DIAM = 6.15;
        VAF.FIN.LEAD_EDGE = 94;
        VAF.FIN.COUNT = 4;
        VAF.FIN.LAMBDA_VEC = [0 90 180 270 360];
        VAF.FIN.C_ROOT = 13.5;
        VAF.FIN.C_TIP = 2.625;
        VAF.FIN.SPAN = 6.5;
        VAF.FIN.SWEEP = 9.625;
        VAF.FIN.THICKNESS = 0.185;
        VAF.FIN.ROUNDED_LEADING_EDGE = 0; % square/rectangular edge chosen
        FIN_CP_CALC('Trapezoid');
        VAF.FIN.AREA_REF = (0.5 * VAF.FIN.BODY_DIAM)^2;
        VAF.FIN.AREA_REF_DRAG = VAF.FIN.COUNT * VAF.FIN.THICKNESS * VAF.FIN.SPAN;
        VAF.FIN.COS_GAMMA = cos(atan((VAF.FIN.SWEEP + 0.5*(VAF.FIN.C_TIP - VAF.FIN.C_ROOT))/VAF.FIN.SPAN));
        FIN_MACH_INTERP();
        VAF.FIN.INTERFERE_CORR = [1, 1.0;...
                                  2, 1.0;...
                                  3, 1.0;...
                                  4, 1.0;...
                                  5, 0.948;...
                                  6, 0.913;...
                                  7, 0.854;...
                                  8, 0.810;...
                                  9, 0.750];
        
        %% Air Properties
        
        VAF.GAS.GAMMA_VS_TEMP = [-40	1.401;...
            -20	1.401;...
            0	1.401;...
            10	1.401;...
            20	1.401;...
            30	1.401;...
            40	1.401;...
            50	1.401;...
            60	1.401;...
            70	1.401;...
            80	1.400;...
            90	1.400;...
            100	1.400;...
            120	1.400;...
            140	1.399;...
            160	1.399;...
            180	1.399;...
            200	1.398;...
            300	1.394;...
            400	1.389;...
            500	1.383;...
            750	1.367;...
            1000	1.351;...
            1500	1.329];
end

% return VAF;

    function NC_CREATE(type, DIAM, LENGTH)
        % creates nosecone variables needed for CP calculation
        % assumes nosecone subsonic pressure drag is zeroed out for now,
        % may need fixing later if time permits
        switch type
            case 'Conical'
                RAD_CHAR = DIAM * 0.5;
                X_VAL = linspace(0, LENGTH, 1000);
                Y_VAL_FUNC = @(X) X*RAD_CHAR / LENGTH;
                
                % DRAG CHARACTERISTICS
                VAF.NC.FINENESS_RATIO = 0.5*DIAM / LENGTH;
                VAF.NC.CD_PRESSURE_SUBSONIC = 0;
                VAF.NC.CD_TRANSONIC = 1.0 * sin(atan(0.5*DIAM / LENGTH));
                VAF.NC.CD_SUPERSONIC = (2.1 / (1 + 4 * (0.5 * DIAM / LENGTH)^2)) + (0.5 / sqrt((1 + 4 * (0.5 * DIAM / LENGTH)^2)*(1.3^2 - 1)));
                VAF.NC.DERIV_CD_MA_TRANSONIC = (4/2.4)*(1-0.5*VAF.NC.CD_TRANSONIC);
                VAF.NC.CD_SUBSONIC_HIGH_PARAM_A = VAF.NC.CD_TRANSONIC - CD_PRESSURE_SUBSONIC;
                VAF.NC.CD_SUBSONIC_HIGH_PARAM_B = VAF.NC.DERIV_CD_MA_TRANSONIC / VAF.NC.CD_SUBSONIC_HIGH_PARAM_A;
                TRANSONIC_MAT = [1.3^2, 1.3, 1;...
                    1, 1, 1;...
                    2, 1, 0];
                TRANSONIC_B = [VAF.NC.CD_SUPERSONIC; VAF.NC.CD_TRANSONIC; VAF.DERIV_CD_MA_TRANSONIC];
                VAF.NC.TRANSONIC_CONSTANTS = TRANSONIC_MAT\TRANSONIC_B;
                
                VAF.NC.TYPE = type;
                VAF.NC.CONE_CORRECTION = 1.0;

            case 'Tangent Ogive'
                RAD_CHAR = ((0.5 * DIAM)^2 + LENGTH^2)/(2 * (DIAM*0.5));
                X_VAL = linspace(0, LENGTH, 1000);
                Y_VAL_FUNC = @(X) (sqrt(RAD_CHAR^2 - (LENGTH - X).^2) + (0.5*DIAM) - RAD_CHAR);
                VAF.NC.CONE_CORRECTION = 0.72 + 0.82;
                
                % DRAG CHARACTERISTICS
                VAF.NC.FINENESS_RATIO = 0.5*DIAM / LENGTH;
                VAF.NC.CD_PRESSURE_SUBSONIC = 0;
                VAF.NC.CD_TRANSONIC = 1.0 * sin(atan(0.5*DIAM / LENGTH));
                VAF.NC.CD_SUPERSONIC = (2.1 / (1 + 4 * (0.5 * DIAM / LENGTH)^2)) + (0.5 / sqrt((1 + 4 * (0.5 * DIAM / LENGTH)^2)*(1.3^2 - 1)));
                VAF.NC.DERIV_CD_MA_TRANSONIC = (4/2.4)*(1-0.5*VAF.NC.CD_TRANSONIC);
                VAF.NC.CD_SUBSONIC_HIGH_PARAM_A = VAF.NC.CD_TRANSONIC - CD_PRESSURE_SUBSONIC;
                VAF.NC.CD_SUBSONIC_HIGH_PARAM_B = VAF.NC.DERIV_CD_MA_TRANSONIC / VAF.NC.CD_SUBSONIC_HIGH_PARAM_A;
                TRANSONIC_MAT = [1.3^2, 1.3, 1;...
                    1, 1, 1;...
                    2, 1, 0];
                TRANSONIC_B = [VAF.NC.CD_SUPERSONIC; VAF.NC.CD_TRANSONIC; VAF.DERIV_CD_MA_TRANSONIC];
                VAF.NC.TRANSONIC_CONSTANTS = INV(TRANSONIC_MAT, TRANSONIC_B);
                
                VAF.NC.TYPE = type;
                VAF.NC.CONE_CORRECTION = 0.72 + 0.82;
                

            case 'Elliptical'
                RAD_CHAR = DIAM * 0.5;
                X_VAL = linspace(0, LENGTH, 1000);
                Y_VAL_FUNC = @(X) RAD_CHAR * sqrt(1 - (1 - (X/LENGTH))^2);
                VAF.NC.CD_PRESSURE_SUBSONIC = 0;
                VAF.NC.TYPE = type;

%             case 'Parabolic'
%                 RAD_CHAR = DIAM * 0.5;
%                 X_VAL = linspace(0, LENGTH, 1000);
%                 Y_VAL_FUNC = @(X) RAD_CHAR * sqrt(1 - (1 - (X/LENGTH))^2);
% 

        end
        
        LEAD_EDGE_AREA = 0;
        TRAIL_EDGE_AREA = pi * Y_VAL_FUNC(X_VAL(end))^2;
        AREA_FUNC = @(X) (Y_VAL_FUNC(X)).^2 * pi;
        VOL = integral(AREA_FUNC, 0, X_VAL(end));
        VAF.NC.X_CP = BODY_CP_CALC(LENGTH, TRAIL_EDGE_AREA, VOL, LEAD_EDGE_AREA);
        VAF.NC.AREA_REF = pi * (DIAM * 0.5)^2;
        VAF.NC.AREA_PLAN = 2 * integral(Y_VAL_FUNC, 0, X_VAL(end));
        VAF.NC.AREA_WET = 2*pi*integral(Y_VAL_FUNC, 0, X_VAL(end));

    end

    function CP = BODY_CP_CALC(LENGTH, AREA_TRAIL_EDGE, VOL, AREA_LEAD_EDGE)
        % calculates CP for body components
        CP = (LENGTH * AREA_TRAIL_EDGE - VOL)/(AREA_TRAIL_EDGE - AREA_LEAD_EDGE);
        if isnan(CP)
            CP = 0;
        end
    end

    function FIN_CP_CALC(type)
        switch type
            case 'Trapezoid'
               VAF.FIN.X_F_SUBSONIC = (VAF.FIN.SWEEP / 3) * ( VAF.FIN.C_ROOT + 2 * VAF.FIN.C_TIP)/(VAF.FIN.C_ROOT + VAF.FIN.C_TIP) + (1/6)*(VAF.FIN.C_ROOT^2 + VAF.FIN.C_TIP^2 + VAF.FIN.C_ROOT*VAF.FIN.C_TIP)/(VAF.FIN.C_ROOT + VAF.FIN.C_TIP);
               VAF.FIN.AREA = 0.5 * (VAF.FIN.C_ROOT + VAF.FIN.C_TIP) * VAF.FIN.SPAN;
               VAF.FIN.MAC = integral(@(x) (((VAF.FIN.C_ROOT - VAF.FIN.C_TIP)/(0 - VAF.FIN.SPAN))*(x)).^2, 0, VAF.FIN.SPAN)/VAF.FIN.AREA;
               VAF.FIN.AREA_WET_DOUBLE_TOTAL = VAF.FIN.AREA * 2 * VAF.FIN.COUNT;
       end
        
    end

    function FIN_MACH_INTERP()
        Mgre2 = @(x) (2*VAF.FIN.SPAN^2 * sqrt(x^2 - 1)/VAF.FIN.AREA - 0.67)/(2 * (2 * VAF.FIN.SPAN^2 * sqrt(x^2 - 1))/VAF.FIN.AREA - 1);
        dx = 0.0005;
        deriv = (Mgre2(2+dx) - Mgre2(2-dx))/(2 * dx);
        VAF.FIN.p_arr = [0.25; 0; Mgre2(2); deriv; 0; 0];
        VAF.FIN.p_mat = [NTH_DERIV(0.5, 5, 0); NTH_DERIV(0.5, 5, 1); NTH_DERIV(2, 5, 0); NTH_DERIV(2, 5, 1); NTH_DERIV(2,5,2); NTH_DERIV(2, 5, 3)];
        VAF.FIN.CP_COEFFS = VAF.FIN.p_mat\VAF.FIN.p_arr;
    end

    function array = NTH_DERIV(val, order, deriv)
       power = (order:-1:0) - deriv;
       coeff = ones(1, order+1);
       for i = 1:deriv
            coeff = coeff .* ([order:-1:0] - (i-1));
       end
       array = ((val * ones(1, order+1)) .^ power).* coeff;
    end
end

