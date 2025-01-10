% Metas.UncLib.Matlab.LinProp V2.7.1
% Michael Wollensack METAS - 24.05.2023
% Dion Timmermann PTB - 22.06.2022
%
% This class supports the creation of uncertainty objects and subsequent
% calculation with them as well as storage of the results. It can handle
% complex-valued and multivariate quantities. Internally, LinProp objects
% in MATLAB are wrappers of the .NET interface of the <a href="https://www.metas.ch/unclib">METAS UncLib library</a>.
%
% Commonly used Constructors  (Round brackets indicate vectors.)
%   u = LinProp(value)
%   u = LinProp(value, standard_unc, [description])
%   u = LinProp(value, (covariance), [description])
%  (u)= LinProp((value), (covariance), [description])
%   u = LinProp((samples), 'samples', [description], [probability])
%   u = LinProp((samples), 'randomchoices', [description], [probability])
%   u = LinProp(distribution, [id], [description])
%   u = LinProp(value, (sys_inputs), (sys_sensitivities), 'system')
%   Documentation of all constructors available with <a href="matlab: help LinProp/LinProp -displayBanner">help LinProp/LinProp</a>.
%
% Properties and Common Functions
%   The values of u can be accessed through u<a href="matlab:help LinProp.Value -displayBanner">.Value</a> or <a href="matlab:help LinProp.get_value -displayBanner">get_value</a>(u), 
%   and the standard uncertainties through u<a href="matlab:help LinProp.StdUnc -displayBanner">.StdUnc</a> or <a href="matlab:help LinProp.get_stdunc -displayBanner">get_stdunc</a>(u).
%   Many common MATLAB functions are available, see <a href="matlab:methods(LinProp(1))">list of all methods</a>.
%
% Uncertainty Methods
%   <a href="matlab:help LinProp.get_correlation -displayBanner"        >get_correlation</a>         Correlation matrix
%   <a href="matlab:help LinProp.get_covariance -displayBanner"         >get_covariance</a>          Covariance matrix
%   <a href="matlab:help LinProp.get_coverage_interval -displayBanner"  >get_coverage_interval</a>   Coverage interval bounds
%   <a href="matlab:help LinProp.get_jacobi -displayBanner"             >get_jacobi</a>              Uncertainty contributions of base inputs
%   <a href="matlab:help LinProp.get_unc_component -displayBanner"      >get_unc_component</a>       Uncertainty contributions of intermediate results
%   <a href="matlab:help LinProp.get_jacobi2 -displayBanner"            >get_jacobi2</a>             Sensitivities to intermediate results
%   <a href="matlab:help LinProp.get_idof -displayBanner"               >get_idof</a>                Inverse degrees of freedom
%   <a href="matlab:help LinProp.get_moment -displayBanner"             >get_moment</a>              n'th central moment
%   <a href="matlab:help LinProp.unc_budget -displayBanner"             >unc_budget</a>              Opens budget window (LinProp only)
%
% Interpolation and Integration Methods
%   <a href="matlab:help LinProp.integrate -displayBanner"          >integrate</a>           Integration with cumulative result
%   <a href="matlab:help LinProp.integrate2 -displayBanner"         >integrate2</a>          Integration with scalar result
%   <a href="matlab:help LinProp.interpolation -displayBanner"      >interpolation</a>       Interpolation
%   <a href="matlab:help LinProp.interpolation2 -displayBanner"     >interpolation2</a>      Interpolation with linear unc. propagation
%   <a href="matlab:help LinProp.spline -displayBanner"             >spline</a>              Spline interpolation
%   <a href="matlab:help LinProp.spline2 -displayBanner"            >spline2</a>             Spline interpolation with linear unc. propagation
%   <a href="matlab:help LinProp.splinecoefs -displayBanner"        >splinecoefs</a>         Coefficients of interpolation spline
%   <a href="matlab:help LinProp.splineintegrate -displayBanner"    >splineintegrate</a>     Spline integration with cumulative result
%   <a href="matlab:help LinProp.splineintegrate2 -displayBanner"   >splineintegrate2</a>    Spline integration with scalar result
%
% Object Behavior
%   Scalar LinProp objects behave like MATLAB fundamental types with
%   respect to copy operations. Copies are independent values. Operations
%   that you perform on one object do not affect copies of that object.
%   Non-scalar LinProp objects are referenced by their handle variable.
%   Copies of the handle variable refer to the same object. Operations that
%   you perform on a handle object are visible from all handle variables
%   that reference that object.
%
%   B = <a href="matlab:help LinProp.copy -displayBanner">copy</a>(A) copies each element in the array of handles A to a new
%   array of handles B.
% 
% See also LinProp/LinProp.


classdef LinProp
    properties
        NetObject
    end
    properties (SetAccess = private)
        Value
        StdUnc
        IsComplex
        IsArray
    end
    methods
        function obj = LinProp(varargin)
            % Constructor to create LinProp uncertainty objects.
            %
            % u = LinProp(value) creates an uncertainty object without any
            % uncertainties but the specified value. value can real- or complex-valued
            % and have any shape. This operation can be used to preallocate variables.
            %
            % u = LinProp(value, standard_unc, [description]) creates a real-valued
            % scalar uncertainty object with the specified value and standard
            % uncertainty. value and standard_unc must be real-valued scalars.
            % Optionally, a description can be specified (see below).
            %
            % u = LinProp(value, covariance, [description]) creates an uncertainty
            % object based on a complex and/or non-scalar value and an associated
            % covariance matrix. value must be a complex scalar or a real- or
            % complex-valued matrix. covariance must be a square, real-valued matrix.
            % The rows and columns of covariance match the elements of value interpred
            % as a vector, i.e. value(:). If value is complex-valued the first row and
            % column of covariance relate to real(value(1)), the second row and colum
            % of covariance to imag(value(1)), the third row and colum of covariance to
            % real(value(2)), etc. Optionally, a description can be specified (see
            % below).
            %
            % u = LinProp(samples, 'samples', [description], [probability]) creates a
            % scalar or vector uncertainty object from a column vector or matrix of
            % real- or complex-valued samples. u will be a scalar or column vector. The
            % number of elements in u is the same as the number of columns of samples.
            % samples must have at least one more row than columns. Optionally, a
            % description can be specified (see below).
            %
            % u = LinProp(value, sys_inputs, sys_sensitivities, 'system') creates an
            % uncertainty object with the specified value and sensitivities to the
            % specified system inputs. This approach can be used to bridge a function,
            % e.g. a numerical method like nonlinear least squares. value must be the
            % result of this function, sys_inputs a vector of uncertainty objects used
            % as inputs to this function, and sys_sensitivities a vector of sensitivies
            % of value to changes of sys_inputs. This constructor is LinProp only!
            % 
            % Use u = LinProp(xml_string), u = LinProp(filepath, 'xml_file'), or 
            % u = LinProp(filepath, 'binary_file') to load an uncertainty object which 
            % has previously been exported with xml_string(u), xml_file(u, filepath), 
            % or binary_file(u, filepath).
            % 
            % LinProp(value, standard_unc, idof, [id, description]) creates a scalar,
            % real-valued uncertainty object based on the specified value, standard
            % uncertainy, and inverse degree of freedom. All three variables must be
            % real-valued scalars. Additionally, an id and description can be specified
            % for the uncertainty.
            %
            % With several constructors it is possible define a description to an
            % uncertainty object. This description must be a char array and will later
            % show in the uncertainty budget. Use '\t' to mark levels of hirachy and
            % thus group uncertainties.
            
            % The assemblies are guaranteed to be loaded through the
            % constant UncHelper property.
            switch nargin
                case 1
                    switch class(varargin{1})
                        case 'LinProp'
                            obj = varargin{1};
                        case 'double'
                            if numel(varargin{1}) == 1
                                if ~isreal(varargin{1})
                                    % ComplexUncNumber
                                    temp = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                                    temp.InitDblReIm(real(varargin{1}), imag(varargin{1}));
                                    obj.NetObject = temp;
                                else
                                    % RealUncNumber
                                    obj.NetObject = Metas.UncLib.LinProp.UncNumber(real(varargin{1}));
                                end
                            else
                                v = LinProp.Double2Array(varargin{1});
                                if ~isreal(varargin{1})
                                    % ComplexUncArray
                                    obj.NetObject = LinProp.UncHelper.ComplexUncNArray(v);
                                else
                                    % RealUncArray
                                    obj.NetObject = LinProp.UncHelper.RealUncNArray(v);
                                end
                            end
                        case 'Metas.UncLib.LinProp.UncNumber'
                            obj.NetObject = varargin{1};
                        case 'Metas.UncLib.Core.Real<Metas*UncLib*LinProp*UncNumber>'
                            obj.NetObject = varargin{1}.Item;
                        case 'Metas.UncLib.Core.Complex<Metas*UncLib*LinProp*UncNumber>'
                            obj.NetObject = varargin{1};
                        case 'Metas.UncLib.Core.Ndims.RealNArray<Metas*UncLib*LinProp*UncNumber>'
                            obj.NetObject = varargin{1};
                        case 'Metas.UncLib.Core.Ndims.ComplexNArray<Metas*UncLib*LinProp*UncNumber>'
                            obj.NetObject = varargin{1};
                        case 'char'
                            obj.NetObject = LinProp.XmlString2LinProp(varargin{1}).NetObject;
                        otherwise
                            if isa(varargin{1}, 'Metas.UncLib.Core.Unc.Distribution')
                                obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, UncInputId(), '');
                            else
                                error('Wrong type of input arguments')
                            end
                    end
                case 2
                    if isa(varargin{1}, 'double') && isa(varargin{2}, 'double')
                        if numel(varargin{1}) == 1
                            if ~isreal(varargin{1})
                                % ComplexUncNumber
                                v = LinProp.Double2ComplexNumber(varargin{1});
                                cv = LinProp.Double2Array(varargin{2});
                                obj.NetObject = LinProp.UncHelper.ComplexUncNumber(v, cv.Matrix, 0);
                            else
                                % RealUncNumber
                                obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, varargin{2});
                            end
                        else
                            v = LinProp.Double2Array(varargin{1});
                            cv = LinProp.Double2Array(varargin{2});
                            if ~isreal(varargin{1})
                                % ComplexUncArray
                                obj.NetObject = LinProp.UncHelper.ComplexUncNArray(v, cv.Matrix, 0);
                            else
                                % RealUncArray
                                obj.NetObject = LinProp.UncHelper.RealUncNArray(v, cv.Matrix, 0);
                            end
                        end
                    elseif isa(varargin{1}, 'char') && isa(varargin{2}, 'char')
                        switch lower(varargin{2})
                            case 'xml_file'
                                obj.NetObject = LinProp.XmlFile2LinProp(varargin{1}).NetObject;
                            case 'binary_file'
                                obj.NetObject = LinProp.BinaryFile2LinProp(varargin{1}).NetObject;
                            otherwise
                                error('Wrong file type')
                        end
                    elseif isa(varargin{1}, 'double') && isa(varargin{2}, 'char')
                        switch lower(varargin{2})
                            case 'samples'
                                s = LinProp.Double2Array(varargin{1});
                                if size(varargin{1}, 2) == 1
                                    if ~isreal(varargin{1})
                                        % ComplexUncNumber
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNumberFromSamples(s.Vector);
                                    else
                                        % RealUncNumber
                                        obj.NetObject = LinProp.UncHelper.RealUncNumberFromSamples(s.Vector);
                                    end
                                else
                                    if ~isreal(varargin{1})
                                        % ComplexUncArray
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNArrayFromSamples(s.Matrix);
                                    else
                                        % RealUncArray
                                        obj.NetObject = LinProp.UncHelper.RealUncNArrayFromSamples(s.Matrix);
                                    end
                                end
                            case 'randomchoices'
                                s = LinProp.Double2Array(varargin{1});
                                if size(varargin{1}, 2) == 1
                                    if ~isreal(varargin{1})
                                        % ComplexUncNumber
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNumberFromRandomChoices(s.Vector);
                                    else
                                        % RealUncNumber
                                        obj.NetObject = LinProp.UncHelper.RealUncNumberFromRandomChoices(s.Vector);
                                    end
                                else
                                    if ~isreal(varargin{1})
                                        % ComplexUncArray
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNArrayFromRandomChoices(s.Matrix);
                                    else
                                        % RealUncArray
                                        obj.NetObject = LinProp.UncHelper.RealUncNArrayFromRandomChoices(s.Matrix);
                                    end
                                end
                            otherwise
                                error('Wrong type of input arguments')
                        end
                    elseif isa(varargin{1}, 'Metas.UncLib.Core.Unc.Distribution') && isa(varargin{2}, 'Metas.UncLib.Core.Unc.InputId')
                        obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, varargin{2}, '');
                    elseif isa(varargin{1}, 'Metas.UncLib.Core.Unc.Distribution') && isa(varargin{2}, 'char')
                        obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, UncInputId(), sprintf(varargin{2}));
                    else
                        error('Wrong type of input arguments')
                    end
                case 3
                    if isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && isa(varargin{3}, 'double')
                        if numel(varargin{1}) == 1
                            if ~isreal(varargin{1})
                                % ComplexUncNumber
                                v = LinProp.Double2ComplexNumber(varargin{1});
                                cv = LinProp.Double2Array(varargin{2});
                                obj.NetObject = LinProp.UncHelper.ComplexUncNumber(v, cv.Matrix, varargin{3});
                            else
                                % RealUncNumber
                                obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, varargin{2}, varargin{3});
                            end
                        else
                            v = LinProp.Double2Array(varargin{1});
                            cv = LinProp.Double2Array(varargin{2});
                            if ~isreal(varargin{1})
                                % ComplexUncArray
                                obj.NetObject = LinProp.UncHelper.ComplexUncNArray(v, cv.Matrix, varargin{3});
                            else
                                % RealUncArray
                                obj.NetObject = LinProp.UncHelper.RealUncNArray(v, cv.Matrix, varargin{3});
                            end
                        end
                    elseif isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && isa(varargin{3}, 'char')
                        if numel(varargin{1}) == 1
                            if ~isreal(varargin{1})
                                % ComplexUncNumber (Description)
                                v = LinProp.Double2ComplexNumber(varargin{1});
                                cv = LinProp.Double2Array(varargin{2});
                                obj.NetObject = LinProp.UncHelper.ComplexUncNumber(v, cv.Matrix, UncInputId(), sprintf(varargin{3}));
                            else
                                % RealUncNumber (Description)
                                obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, varargin{2}, 0, UncInputId(), sprintf(varargin{3}));
                            end
                        else
                            v = LinProp.Double2Array(varargin{1});
                            cv = LinProp.Double2Array(varargin{2});
                            if ~isreal(varargin{1})
                                % ComplexUncArray (Description)
                                obj.NetObject = LinProp.UncHelper.ComplexUncNArray(v, cv.Matrix, UncInputId(), sprintf(varargin{3}));
                            else
                                % RealUncArray (Description)
                                obj.NetObject = LinProp.UncHelper.RealUncNArray(v, cv.Matrix, UncInputId(), sprintf(varargin{3}));
                            end
                        end
                    elseif isa(varargin{1}, 'double') && isa(varargin{2}, 'char') && isa(varargin{3}, 'char')
                        switch lower(varargin{2})
                            case 'samples'
                                s = LinProp.Double2Array(varargin{1});
                                if size(varargin{1}, 2) == 1
                                    if ~isreal(varargin{1})
                                        % ComplexUncNumber
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNumberFromSamples(s.Vector, UncInputId(), sprintf(varargin{3}));
                                    else
                                        % RealUncNumber
                                        obj.NetObject = LinProp.UncHelper.RealUncNumberFromSamples(s.Vector, UncInputId(), sprintf(varargin{3}));
                                    end
                                else
                                    if ~isreal(varargin{1})
                                        % ComplexUncArray
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNArrayFromSamples(s.Matrix, UncInputId(), sprintf(varargin{3}));
                                    else
                                        % RealUncArray
                                        obj.NetObject = LinProp.UncHelper.RealUncNArrayFromSamples(s.Matrix, UncInputId(), sprintf(varargin{3}));
                                    end
                                end
                            case 'randomchoices'
                                s = LinProp.Double2Array(varargin{1});
                                if size(varargin{1}, 2) == 1
                                    if ~isreal(varargin{1})
                                        % ComplexUncNumber
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNumberFromRandomChoices(s.Vector, UncInputId(), sprintf(varargin{3}));
                                    else
                                        % RealUncNumber
                                        obj.NetObject = LinProp.UncHelper.RealUncNumberFromRandomChoices(s.Vector, UncInputId(), sprintf(varargin{3}));
                                    end
                                else
                                    if ~isreal(varargin{1})
                                        % ComplexUncArray
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNArrayFromRandomChoices(s.Matrix, UncInputId(), sprintf(varargin{3}));
                                    else
                                        % RealUncArray
                                        obj.NetObject = LinProp.UncHelper.RealUncNArrayFromRandomChoices(s.Matrix, UncInputId(), sprintf(varargin{3}));
                                    end
                                end
                            otherwise
                                error('Wrong type of input arguments')
                        end
                    elseif isa(varargin{1}, 'Metas.UncLib.Core.Unc.Distribution') && isa(varargin{2}, 'Metas.UncLib.Core.Unc.InputId') && isa(varargin{3}, 'char')
                        obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, varargin{2}, sprintf(varargin{3}));
                    else
                        error('Wrong type of input arguments')
                    end
                case 4
                    if isa(varargin{1}, 'double') && isa(varargin{2}, 'LinProp') && isa(varargin{3}, 'double') && isa(varargin{4}, 'char')
                        switch lower(varargin{4})
                            case 'system'
                                if ~isreal(varargin{1}) || ~isreal(varargin{2}) || ~isreal(varargin{3})
                                    error('Value, system inputs, and system sensitivities must be real-valued.');
                                end
                                obj.NetObject = LinProp.System2LinProp(varargin{1}, varargin{2}, varargin{3}).NetObject;
                            otherwise
                                error('Wrong type of input arguments')
                        end
                    elseif isa(varargin{1}, 'double') && isa(varargin{2}, 'char') && isa(varargin{3}, 'char') && isa(varargin{4}, 'double')
                        switch lower(varargin{2})
                            case 'samples'
                                s = LinProp.Double2Array(varargin{1});
                                if size(varargin{1}, 2) == 1
                                    if ~isreal(varargin{1})
                                        % ComplexUncNumber
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNumberFromSamples(s.Vector, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    else
                                        % RealUncNumber
                                        obj.NetObject = LinProp.UncHelper.RealUncNumberFromSamples(s.Vector, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    end
                                else
                                    if ~isreal(varargin{1})
                                        % ComplexUncArray
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNArrayFromSamples(s.Matrix, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    else
                                        % RealUncArray
                                        obj.NetObject = LinProp.UncHelper.RealUncNArrayFromSamples(s.Matrix, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    end
                                end
                            case 'randomchoices'
                                s = LinProp.Double2Array(varargin{1});
                                if size(varargin{1}, 2) == 1
                                    if ~isreal(varargin{1})
                                        % ComplexUncNumber
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNumberFromRandomChoices(s.Vector, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    else
                                        % RealUncNumber
                                        obj.NetObject = LinProp.UncHelper.RealUncNumberFromRandomChoices(s.Vector, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    end
                                else
                                    if ~isreal(varargin{1})
                                        % ComplexUncArray
                                        obj.NetObject = LinProp.UncHelper.ComplexUncNArrayFromRandomChoices(s.Matrix, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    else
                                        % RealUncArray
                                        obj.NetObject = LinProp.UncHelper.RealUncNArrayFromRandomChoices(s.Matrix, UncInputId(), sprintf(varargin{3}), varargin{4});
                                    end
                                end
                            otherwise
                                error('Wrong type of input arguments')
                        end
                    else
                        error('Wrong type of input arguments')
                    end
                case 5
                    if isa(varargin{1}, 'double') && isa(varargin{2}, 'double') && isa(varargin{3}, 'double') && isa(varargin{5}, 'char')
                        if ~isa(varargin{4}, 'Metas.UncLib.Core.Unc.InputId')
                            varargin{4} = UncInputId(varargin{4});
                        end
                        if numel(varargin{1}) == 1
                            if ~isreal(varargin{1})
                                % ComplexUncNumber
                                v = LinProp.Double2ComplexNumber(varargin{1});
                                cv = LinProp.Double2Array(varargin{2});
                                obj.NetObject = LinProp.UncHelper.ComplexUncNumber(v, cv.Matrix, varargin{3}, varargin{4}, sprintf(varargin{5}));
                            else
                                % RealUncNumber
                                obj.NetObject = Metas.UncLib.LinProp.UncNumber(varargin{1}, varargin{2}, varargin{3}, varargin{4}, sprintf(varargin{5}));
                            end
                        else
                            v = LinProp.Double2Array(varargin{1});
                            cv = LinProp.Double2Array(varargin{2});
                            if ~isreal(varargin{1})
                                % ComplexUncArray
                                obj.NetObject = LinProp.UncHelper.ComplexUncNArray(v, cv.Matrix, varargin{3}, varargin{4}, sprintf(varargin{5}));
                            else
                                % RealUncArray
                                obj.NetObject = LinProp.UncHelper.RealUncNArray(v, cv.Matrix, varargin{3}, varargin{4}, sprintf(varargin{5}));
                            end
                        end
                    else
                        error('Wrong type of input arguments')
                    end
                otherwise
                    error('Wrong number of input arguments')
            end
            % Ensure arrays are internally always stored as matrices.
            if LinProp.IsArrayNet(obj.NetObject)
                if obj.NetObject.ndims == 1
                    obj.NetObject.Reshape(int32([obj.NetObject.numel 1]));
                end
            end 
        end
        function str = string(obj)
            
            pm = sprintf(' \xB1 '); % The plus/minus sign coded as unicode 
            % number so this source code file is not dependent on the encoding.
            
            % Using evalc(disp(x)) prints using the current format setting.
            edisp = @(x) strtrim(evalc('disp(x)'));
            
            str = cell(size(obj));
            
            val_real = get_value(real(obj));
            unc_real = get_stdunc(real(obj));
            sign_real = repmat(' ', size(obj));
            sign_real(val_real < 0) = '-';
            val_real = abs(val_real);
            
            if ~obj.IsComplex
                for ii = 1:numel(obj)
                    str{ii} = [sign_real(ii) '(' edisp(val_real(ii)) pm edisp(unc_real(ii)) ')'];
                end
            else          
                val_imag = get_value(imag(obj));
                unc_imag = get_stdunc(imag(obj));
                sign_imag = repmat('+', size(obj));
                sign_imag(val_imag < 0) = '-';
                val_imag = abs(val_imag);
                
                for ii = 1:numel(obj)
                    str{ii} = [sign_real(ii)  '(' edisp(val_real(ii)) pm edisp(unc_real(ii)) ') ' ...
                               sign_imag(ii) ' (' edisp(val_imag(ii)) pm edisp(unc_imag(ii)) ')i'];
                end
            end
            
            % Strings and the string() function were introduced in Matalb
            % 2016b (version 9.1). Return the cellstr for older versions.
            if ~verLessThan('matlab', '9.1')
                str = string(str);
            end
            
        end
        function display(obj)
            name = inputname(1);
            ds = get(0, 'FormatSpacing');
            if isempty(obj)
                if isequal(ds, 'compact')
                    fprintf('%s =\n     []\n', name);
                else
                    fprintf('\n%s =\n\n     []\n\n', name);
                end
            elseif obj.IsArray
                value = get_value(obj);
                unc = get_stdunc(obj);
                dispAsPages([name '.Value'], value, isequal(ds, 'loose'));
                dispAsPages([name '.StdUnc'], unc, isequal(ds, 'loose'));
            else
                if isequal(ds, 'compact')
                    fprintf('%s =\n  %s\n', name, char(string(obj)));
                else
                    fprintf('\n%s =\n\n  %s\n\n', name, char(string(obj)));
                end
            end
        end
        function o = copy(obj)
            if obj.IsArray
                o = LinProp(Copy(obj.NetObject));
            else
                o = obj;
            end
        end
        function index = end(obj, position, numindices)
            if (numindices == 1)
                index = numel(obj);
            else
                index = size(obj, position);
            end
        end
        function l = length(obj)
            if obj.IsArray
                s = double(obj.NetObject.size);
            else
                s = [1 1];
            end
            l = max(s);
        end
        function n = ndims(obj)
            if obj.IsArray
                n = max(2, double(obj.NetObject.ndims));
            else
                n = 2;
            end
        end
        function n = numel(obj)
            if obj.IsArray
                n = double(obj.NetObject.numel);
            else
                n = 1;
            end
        end
        function e = isempty(obj)
            if obj.IsArray && obj.NetObject.numel == 0
                e = true;
            else
                e = false;
            end
        end
        function varargout = size(obj, varargin)
            %SIZE   Size of array.  
            %   D = SIZE(X), for M-by-N matrix X, returns the two-element row vector
            %   D = [M,N] containing the number of rows and columns in the matrix.
            %   For N-D arrays, SIZE(X) returns a 1-by-N vector of dimension lengths.
            %   Trailing singleton dimensions are ignored.
            %
            %   [M,N] = SIZE(X) for matrix X, returns the number of rows and columns in
            %   X as separate output variables. 
            %   
            %   [M1,M2,M3,...,MN] = SIZE(X) for N>1 returns the sizes of the first N 
            %   dimensions of the array X.  If the number of output arguments N does
            %   not equal NDIMS(X), then for:
            %
            %   N > NDIMS(X), SIZE returns ones in the "extra" variables, i.e., outputs
            %                 NDIMS(X)+1 through N.
            %   N < NDIMS(X), MN contains the product of the sizes of dimensions N
            %                 through NDIMS(X).
            %
            %   M = SIZE(X,DIM) returns the lengths of the specified dimensions in a 
            %   row vector. DIM can be a scalar or vector of dimensions. For example, 
            %   SIZE(X,1) returns the number of rows of X and SIZE(X,[1 2]) returns a 
            %   row vector containing the number of rows and columns.
            %
            %   M = SIZE(X,DIM1,DIM2,...,DIMN) returns the lengths of the dimensions
            %   DIM1,...,DIMN as a row vector.
            %
            %   [M1,M2,...,MN] = SIZE(X,DIM) OR [M1,M2,...,MN] = SIZE(X,DIM1,...,DIMN)
            %   returns the lengths of the specified dimensions as separate outputs.
            %   The number of outputs must equal the number of dimensions provided.
            %
            
            % Write size of all dimensions to s.
            if LinProp.IsArrayNet(obj.NetObject)
                netSize = obj.NetObject.size; % Using a temp variable saves a lot of time.
                s = double(netSize);
            else
                s = [1 1];
            end
            
            % Write all requested dimensions to dims
            if nargin == 1
                % Special case for nargout ~= length(s) if no dims were specificed 
                if nargout > 1
                    if nargout > length(s)
                        s(end+1:nargout) = 1;
                    elseif nargout < length(s)
                        s = [s(1:nargout-1) prod(s(nargout:end))];
                    end
                end
            elseif nargin == 2
                dims = varargin{1};

                % Check if requested dims are valid
                if any(dims < 1 | ceil(dims) ~= dims | isinf(dims))
                    error('Dimension argument must be a positive integer scalar or a vector of positive integers.'); 
                end

                % Add singleton dimensions and reduce s to selected dims
                s = [s ones(1, max(dims)-length(s))];
                s = s(dims);
            else
                if any(cellfun(@(x) ~isscalar(x) || ~isnumeric(x), varargin))
                    error('Dimension argument must be a positive integer scalar within indexing range.');
                end
                dims = cell2mat(varargin);

                % Check if requested dims are valid
                if any(dims < 1 | ceil(dims) ~= dims | isinf(dims))
                    error('Dimension argument must be a positive integer scalar or a vector of positive integers.'); 
                end

                % Add singleton dimensions and reduce s to selected dims
                s = [s ones(1, max(dims)-length(s))];
                s = s(dims);
            end
            
            if nargout == 0 || nargout == 1
                varargout = {s};
            elseif nargout == numel(s)
                varargout = num2cell(s);
            else
                error('Incorrect number of output arguments. Number of output arguments must equal the number of input dimension arguments.');
            end
            
        end
        function y = reshape(x, varargin)
            %RESHAPE Reshape array.
            %   RESHAPE(X,M,N) or RESHAPE(X,[M,N]) returns the M-by-N matrix
            %   whose elements are taken columnwise from X. An error results
            %   if X does not have M*N elements.
            %
            %   RESHAPE(X,M,N,P,...) or RESHAPE(X,[M,N,P,...]) returns an
            %   N-D array with the same elements as X but reshaped to have
            %   the size M-by-N-by-P-by-.... The product of the specified
            %   dimensions, M*N*P*..., must be the same as NUMEL(X).
            %
            %   RESHAPE(X,...,[],...) calculates the length of the dimension
            %   represented by [], such that the product of the dimensions
            %   equals NUMEL(X). The value of NUMEL(X) must be evenly divisible
            %   by the product of the specified dimensions. You can use only one
            %   occurrence of [].
            %
            unknownDimensions = [];
            if numel(varargin) > 1
                unknownDimensions = cellfun(@isempty, varargin);
                if sum(unknownDimensions) > 1
                    error('Size can only have one unknown dimension.');
                elseif sum(unknownDimensions) == 1
                    % Temporarily replace empty dimension with 1, so
                    %   prod() can be used, as the unknown dimension is for now assumed to be 1,
                    %   numel(varargin) == numel(cell2mat(varargin)), and
                    %   checks on numeric values do not fail.
                    varargin(unknownDimensions) = {1};
                    % Correct value will be calculated after other arguments habe been checked.
                end
                if any(not(cellfun(@isscalar, varargin)))
                    error('Size arguments must be integer scalars.');
                end
                s = cell2mat(varargin);
            else
                s = double(varargin{1});
            end
            if numel(s) < 2
                error('MATLAB:getReshapeDims:sizeVector', 'Size vector must have at least two elements.');
            end
            if any(not(isreal(s)))
                error('Size argument cannot be complex.');
            end
            % check if size arguments have integer values (rounding has no effect and not inf, nan also fails this test).
            if any(ceil(s)~=s | isinf(s))
                error('Size arguments must be real integers.');
            end
            % Fix and check dimensions
            if sum(unknownDimensions) == 1
                if mod(numel(x), prod(s)) ~= 0
                    error('Product of known dimensions, %i, not divisible into total number of elements, %i.', prod(s), numel(x));
                else
                    s(unknownDimensions) = numel(x) / prod(s);
                end
            else
                if prod(s) ~= numel(x)
                    error('MATLAB:getReshapeDims:notSameNumel', 'Number of elements must not change. Use [] as one of the size inputs to automatically calculate the appropriate size for that dimension.');
                end
            end
            y = copy(x);
            ym = LinProp.Convert2UncArray(y);
            ym.Reshape(int32(s(:)));
            y = LinProp.Convert2LinProp(ym);
        end
        function C = subsasgn(A, S, B)
            %SUBSASGN Subscripted assignment.
            %   A(I) = B assigns the values of B into the elements of A specified by
            %   the subscript vector I.  B must have the same number of elements as I
            %   or be a scalar. 
            %
            %   A(I,J) = B assigns the values of B into the elements of the rectangular
            %   submatrix of A specified by the subscript vectors I and J.  A colon used as
            %   a subscript, as in A(I,:) = B, indicates all columns of those rows
            %   indicated by vector I. Similarly, A(:,J) = B means all rows of columns J.
            %
            %   A(I,J,K,...) = B assigns the values of B to the submatrix of A specified
            %   by the subscript vectors I, J, K, etc. A colon used as a subscript, as in
            %   A(I,:,K) = B, indicates the entire dimension. 
            %
            %   For both A(I,J) = B and the more general multi-dimensional 
            %   A(I,J,K,...) = B, B must be LENGTH(I)-by-LENGTH(J)-by-LENGTH(K)-... , or
            %   be shiftable to that size by adding or removing singleton dimensions, or
            %   contain a scalar, in which case its value is replicated to form a matrix
            %   of that size.
        
            if any(strcmp('.', {S.type}))
                error('Dot indexing is not supported for variables of this type.');
            elseif any(strcmp('{}', {S.type}))
                error('Brace indexing is not supported for variables of this type.');
            elseif length(S) > 1
                error('Invalid array indexing.');    % This type of error should never appear.
            end
            
            % I describes the index-region of A that values are assigned
            % to. I might be larger than A. In that case, A is extended.
            I = S.subs;
            dimI = numel(I);
            
            % Convert logical indexes to subscripts and check subscripts
            for ii = dimI:-1:1
                if islogical(I{ii})
                    I{ii} = find(I{ii});
                else
                    v = I{ii}(:);
                    if any(ceil(v)~=v | isinf(v) | v <= 0)
                        error('Array indices must be positive integers or logical values.');
                    end
                end
            end
            
            newA = strcmp(class(A), 'double'); %#ok<STISA>
            
            sizeA = size(A);
            numelA = prod(sizeA);
            sizeB = size(B);
            numelB = prod(sizeB);
            isemptyB = (numelB == 0);
            isscalarB = (numelB == 1);
            
            % Special case of null assignment to remove elements
            if isemptyB && isa(B, 'double')
                if sum(~strcmp(I, ':')) > 1
                    error('A null assignment can have only one non-colon index.');
                else
                    if dimI == 1
                        if strcmp(I, ':')
                            C = [];
                            return;
                        else
                            S.subs{1} = true(size(A));
                            S.subs{1}(I{1})= false;
                            if sum(sizeA > 1) == 1 % Is vector for arbitrary ndims
                                C = subsref(A, S);
                            else
                                C = subsref(A, S)';
                            end
                            return;
                        end
                    else
                        dim = find(~strcmp(I, ':'));
                        S.subs{dim} = true(1, size(A, dim));
                        S.subs{dim}(I{dim}) = false;
                        C = subsref(A, S);
                        return;
                    end
                end
            end
            
            % Typecasts
            if ~isa(A, 'LinProp')
                A = LinProp(A);
            end
            if ~isa(B, 'LinProp')
                B = LinProp(B);
            end
            isComplexA = A.IsComplex;
            isComplexB = B.IsComplex;
            isComplex = isComplexA || isComplexB;
            if isComplexA && ~isComplexB
                B = complex(B);
            elseif ~isComplexA && isComplexB
                A = complex(A);
            end
            
            % Replace ':' placeholders 
            % Note: The last dimension can always be used to address
            % all following dimensions.
            if all(sizeA == 0)
                % If A has not been defined yet, the dots (:) refer to the
                % size of B.
                if numel(sizeB) ~= sum(cellfun(@numel, I)>1 | strcmp(I, ':'))    % Singleton dimensions of B are ignored, except the dimensions already match.
                    sizeB_reduced = sizeB(sizeB>1);
                    sizeB_reduced = [sizeB_reduced ones(1, numel(I)-numel(sizeB_reduced))];
                else
                    sizeB_reduced = sizeB;
                end
                tmpProd = 1;
                idx = 1;
                if any(strcmp(I, ':'))
                    if dimI < sum(sizeB_reduced>1)
                        error('Unable to perform assignment because the indices on the left side are not compatible with the size of the right side.');
                    end
                    for ii = 1:(dimI-1)  % Dimensions except the last one
                        if strcmp(I{ii}, ':')
                            I{ii} = 1:sizeB_reduced(idx);
                            tmpProd = tmpProd * sizeB_reduced(idx);
                            idx = idx + 1;
                        elseif numel(I{ii}) > 1
                            if numel(I{ii}) ~= sizeB_reduced(idx)
                                error('Unable to perform assignment because the indices on the left side are not compatible with the size of the right side.');
                            end
                            tmpProd = tmpProd * sizeB_reduced(idx);
                            idx = idx + 1;
                        end
                    end
                    if strcmp(I{dimI}, ':') % Special case for last dimension
                        I{dimI} = 1:(numelB/tmpProd);
                    elseif numel(I{dimI}) ~= numelB/tmpProd
                        error('Unable to perform assignment because the indices on the left side are not compatible with the size of the right side.');
                    end
                end
            else
                % If A has already been defined, the dots (:) refer to the
                % size of A.
                for ii = 1:(dimI-1)  % Dimensions except the last one
                    if strcmp(I{ii}, ':')
                        I{ii} = 1:sizeA(ii);
                    end
                end
                if strcmp(I{dimI}, ':') % Special case for last dimension
                    I{dimI} = 1:(numelA/prod(sizeA(1 : (dimI-1))));   
                end
            end
            
            for ii = numel(I):-1:1
                I_isempty = isempty(I{ii});
                if I_isempty
                    I_maxIndex(ii) = 0;
                else
                    I_maxIndex(ii) = double(max(I{ii}));
                end
            end

            % Assignment of no elements to an empty/new object.
            if any(I_isempty) && numelA == 0
                if numelB <= 1
                    s = I_maxIndex;
                    if numel(s) < 2
                        if isemptyB && newA
                            s = [ones(1,2-numel(s)) s];
                        else
                            s = [s zeros(1,2-numel(s))];
                        end
                    else
                        lastNonSingletonDimension = find(s~=1, 1, 'last');
                        s = s(1:max(2, lastNonSingletonDimension));
                    end
                    C = reshape(LinProp([]), s);
                    return;
                else 
                    error('Unable to perform assignment because the indices on the left side are not compatible with the size of the right side.');
                end
            
            % Linear indexing
            elseif dimI == 1
                % Linear indexing follows some specific rules
                
                if ~isscalarB && numel(I{1}) ~= numelB
                    error('Unable to perform assignment because the left and right sides have a different number of elements.');
                end
                
                % Grow vector if necessary
                if I_maxIndex > numelA
                    if numelA == 0
                        A = zeros(1, I_maxIndex, 'LinProp');
                        if isComplex
                            A = complex(A);
                        end
                    elseif isrow(A)
                        A = [A, LinProp(zeros(1, I_maxIndex-numelA))];
                    elseif iscolumn(A)
                        A = [A; LinProp(zeros(I_maxIndex-numelA, 1))];
                    else
                        error('Attempt to grow array along ambiguous dimension.');
                    end
                end
                
                % Call core library functions to copy values
                am = LinProp.Convert2UncArray(A);
                bm = LinProp.Convert2UncArray(B);
                dest_index = LinProp.IndexMatrix(I);

                if isscalarB
                    am.SetSameItem1d(int32(dest_index - 1), bm.GetItem1d(0));
                else
                    am.SetItems1d(int32(dest_index - 1), bm.GetItems1d(int32(0 : numelB-1)));
                end
                
                C = LinProp.Convert2LinProp(am);
                return;
                
            % Or subscript indexing / partial linear indexing
            else
  
                if dimI < numel(sizeA)
                    % partial linear indexing
                    if max(I{end}) > prod(sizeA(dimI:end))
                        error('Attempt to grow array along ambiguous dimension.');
                    end
                else
                    % Ignore empty and singleton dimensions that have been
                    % indexed but do not exist anyways.
                    for ii = dimI:-1:1
                        I_issingleton(ii) = all(I{ii}(:) == 1);
                    end
                    I_lastRelevant = [find(not(I_isempty | I_issingleton), 1, 'last') 2];
                    I_lastRelevant = I_lastRelevant(1);
                    I = I(1:min(dimI, max(numel(sizeA), I_lastRelevant)));
                    dimI = numel(I);
                    I_maxIndex = I_maxIndex(1:dimI);
                end
                
                % Check dimensions
                if ~isscalarB
                    for ii = dimI:-1:1
                        I_numel(ii) = numel(I{ii});
                    end
                    
                    sizeI_reduced = I_numel(I_numel ~= 1);
                    sizeB_reduced = sizeB(sizeB ~= 1);
                    if ~isequal(sizeI_reduced, sizeB_reduced) && not(any(I_numel == 0) && any(sizeB == 0))
                        error('Unable to perform assignment because the size of the left side is %s and the size of the right side is %s.', ...
                        strjoin(string(I_numel), '-by-'), ...
                        strjoin(string(sizeB), '-by-'));
                    end
                    
                end
                    
                % Expand A, if the addressed area is larger
                if numel(I_maxIndex) > numel(sizeA)
                    sizeA(end+1:numel(I_maxIndex)) = 0; % Expand size vector for A, if nI is larger
                end
                sA_nI = [sizeA(1 : (dimI-1)), prod(sizeA(dimI:end))]; % size of A, when using the same number of dimensions as nI;
                if any(I_maxIndex > sA_nI)
                    A2 = LinProp(zeros(max(I_maxIndex, sA_nI)));
                    if isComplex
                        A2 = complex(A2);
                    end
                    if numel(A) == 0
                        A = A2;
                    else
                        % Copy over existing values from A to A2...
                        A = subsasgn(A2, substruct('()', arrayfun(@(x) (1:x), size(A), 'UniformOutput', false)), A);
                    end
                end
                
                % Remove trailing singleton dimensions that might have been
                % addressed. These might actually not exist if the
                % subscript used was 1.
                while numel(I) > 2 && numel(I{end}) == 1 && I{end} == 1
                    I(end) = [];
                end
                
                % Call core library functions to copy values
                am = LinProp.Convert2UncArray(A);
                bm = LinProp.Convert2UncArray(B);
                dest_index = LinProp.IndexMatrix(I);

                if isscalarB
                    am.SetSameItemNd(int32(dest_index - 1), bm.GetItem1d(0));
                else
                    src_subs = arrayfun(@(x) 1:x, sizeB, 'UniformOutput', false);
                    src_index  = LinProp.IndexMatrix(src_subs);

                    am.SetItemsNd(int32(dest_index - 1), bm.GetItemsNd(int32(src_index - 1)));
                end

                C = LinProp.Convert2LinProp(am);
                return;

            end
               
        end
        function n = numArgumentsFromSubscript(~, ~, ~)
            % Number of arguments returned by subsref and required by subsasgn. 
            %
            % When addressing a = {1, 2, 3} with a{:}, this function would
            % return 3. When addressing a = {1, 2, 3} with a(:), this
            % function would return 1. 
            %
            % This class does not support brace indexing. When using dot
            % indexing on a LinProp matrix, e.g. a = LinProp(1:3); a.Value, 
            % we want to return one matrix containing all the values of a.
            % Thus, this function always returns 1.
            %
            n = 1;
        end
        function B = subsref(A, S)
            %SUBSREF Subscripted reference.
            %   A(I) is an array formed from the elements of A specified by the
            %   subscript vector I.  The resulting array is the same size as I except
            %   for the special case where A and I are both vectors.  In this case,
            %   A(I) has the same number of elements as I but has the orientation of A.
            %
            %   A(I,J) is an array formed from the elements of the rectangular
            %   submatrix of A specified by the subscript vectors I and J.  The
            %   resulting array has LENGTH(I) rows and LENGTH(J) columns.  A colon used
            %   as a subscript, as in A(I,:), indicates all columns of those rows
            %   indicated by vector I. Similarly, A(:,J) = B means all rows of columns
            %   J.
            %
            %   For multi-dimensional arrays, A(I,J,K,...) is the subarray specified by
            %   the subscripts.  The result is LENGTH(I)-by-LENGTH(J)-by-LENGTH(K)-...
            
            if strcmp('.', S(1).type)
                B = builtin('subsref', A, S);
            elseif strcmp('{}', S(1).type)
                error('Brace indexing is not supported for variables of this type.');
            else

                ni = numel(S(1).subs);
                if ni == 0
                    B = copy(A);
                else

                    sizeA = size(A);
                    isvectorA = sum(sizeA > 1) == 1;
                    src_subs = S(1).subs;
                    sizeB = [];

                    % Convert logical indexes to subscripts
                    for ii = 1:ni
                        if islogical(src_subs{ii})
                            src_subs{ii} = find(src_subs{ii});
                        end
                    end

                    % This is a very special case. If linear indexing is used, but
                    % the linear indexes are arranged in form of a matrix, the
                    % output has the shape of the matrix. This does not apply to
                    % logical indexes.
                    if ni == 1 && ~isvector(src_subs{1})
                        sizeB = int32(size(src_subs{1}));   % Save shape of output for later.
                        src_subs{1} = src_subs{1}(:);       % But conform to vector for processing.
                    end

                    sizeA_extended = [sizeA ones(1, ni-numel(sizeA))];
                    % Replace ':' placeholders and ensure indexes are integers.
                    % Note: The last dimension can always be used to address
                    % all following dimensions.
                    for ii = 1:(ni-1)  % Dimensions except the last one
                        if strcmp(src_subs{ii}, ':')
                            src_subs{ii} = 1:sizeA_extended(ii);
                        else
                            originalValue = src_subs{ii};
                            src_subs{ii} = int32(src_subs{ii});
                            if ~isequal(originalValue, double(src_subs{ii}))
                                error('Array indices must be positive integers or logical values.');
                            end
                        end
                    end
                    if strcmp(src_subs{ni}, ':') % Special case for last dimension
                        src_subs{ni} = (1:(numel(A)/prod(sizeA_extended(1 : (ni-1)))))';
                    else
                        originalValue = src_subs{ni};
                        src_subs{ni} = int32(src_subs{ni});
                        if ~isequal(originalValue, double(src_subs{ni}))
                            error('Array indices must be positive integers or logical values.');
                        end
                    end

                    % Handling (partial) linear indexing
                    if ni == 1 && isvectorA
                        % If A is a vector and indexed by a vector, the output has the same shape as A. 
                        % This does not apply if the index is ':'.
                        if ~strcmp(S(1).subs{1}, ':') && isempty(sizeB)
                            sizeB = int32(sizeA);
                            sizeB(sizeA > 1) = int32(numel(src_subs{1}));
                        end
                    else
                        % Determine the size of A based on the used
                        % indexing and reshape if necessary.
                        sizeAnew = [sizeA_extended(1:ni-1) prod(sizeA_extended(ni:end))];
                        if numel(sizeAnew) == 1
                            % If linear indexing is used, the shape of the
                            % output is determined by the shape of the index.
                            if isrow(src_subs{1}) 
                                sizeAnew = [1 sizeAnew(1)];
                                sizeB    = int32([1 numel(src_subs{1})]);
                            else % src_subs{1} is a column vector or a matrix(!).
                                sizeAnew = [sizeAnew(1) 1];
                            end
                        end
                        if ~isequal(sizeAnew, sizeA)
                            A = reshape(A, sizeAnew);
                            sizeA = sizeAnew;
                            isvectorA = sum(sizeA > 1) == 1;
                        end
                    end

                    % If the size of B is not determined by some special
                    % case above, calculate it now
                    if isempty(sizeB)
                        for ii = ni:-1:1
                            sizeB(ii) = int32(numel(src_subs{ii}));
                        end
                        % Trailing singleton dimensions are removed
                        if numel(sizeB) > 2
                            lastNonSingletonDimension = find(sizeB~=1, 1, 'last');
                            if lastNonSingletonDimension < 2
                                sizeB = sizeB(1:2);
                            elseif ~isempty(lastNonSingletonDimension)
                                sizeB = sizeB(1:lastNonSingletonDimension);
                            end
                        end
                    end
                    for ii = numel(sizeB):-1:1
                        dest_subs{ii} = 1:sizeB(ii);
                    end

                    % Create the UncArrays and index matricies for copying
                    am = LinProp.Convert2UncArray(A);
                    src_index  = LinProp.IndexMatrix(src_subs);
                    dest_index = LinProp.IndexMatrix(dest_subs);
                    if A.IsComplex
                       bm = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                       bm.InitNd(sizeB);
                    else
                       bm = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                       bm.InitNd(sizeB);
                    end
                    % If A is a scalar, the UncArray am will have at most 2
                    % dimensions. If A was addressed with more than 2
                    % dimensions, e.g. A(1, 1, 1), we simply ignore the
                    % other dimensions. If the indices were anything other
                    % than 1, A would have been reshaped above to not be a
                    % scalar.
                    if prod(sizeA) == 1 && ni > 2
                        src_index = src_index(:, 1:2);
                    end
                    
                    % Copy the selected elements
                    try
                        if prod(sizeB) == 1
                            if ni == 1
                                B = LinProp(am.GetItem1d(int32(src_index - 1)));
                            else
                                B = LinProp(am.GetItemNd(int32(src_index - 1)));
                            end
                        else
                            % If we reach this point, A is guaranteed to be a
                            % matrix.
                            bm.SetItemsNd(int32(dest_index - 1), am.GetItemsNd(int32(src_index - 1)));
                            B = LinProp(bm);
                        end
                    catch e
                        
                        % Some index was incorrect. Test the subscripts to print typical matlab error messages.
                        if any(cellfun(@(v) any(isinf(v) | v <= 0), src_subs))
                            error('Array indices must be positive integers or logical values.');
                        end
                        if ni == 1 && isvectorA
                            if any(src_subs{1} > numel(A))
                                error('Index exceeds the number of array elements (%i).', numel(A));
                            end
                        else
                            too_large = arrayfun(@(m, v) any(v{1} > m), sizeA(1:ni), src_subs);
                            if any(too_large)
                                error('Index in position %i exceeds array bounds (must not exceed %i).', find(too_large>0, 1), sizeA(find(too_large>0, 1)));
                            end
                        end
                        
                        % Oterhwise rethrow prior error (this should not happen).
                        rethrow(e);
                    end
                    
                end
                
                % after S(1).type == '()' has been processed
                if length(S) > 1
                    B = subsref(B, S(2:end));
                end
            end
        end
        function c = cat(dim, a, varargin)
            
            c = a;
                
            if numel(varargin) > 0
                ndimsA = ndims(a);
                sizeA = size(a);
                sizeVararginElements = cell(1, numel(varargin));
                for ii = 1:numel(varargin)
                    if ndims(varargin{ii}) ~= ndimsA
                        error('Dimensions of arrays being concatenated are not consistent.');
                    end
                    sizeVararginElements{ii} = size(varargin{ii});
                    for kk = 1:ndimsA
                        if kk ~= dim && sizeA(kk) ~= sizeVararginElements{ii}(kk)
                            error('Dimensions of arrays being concatenated are not consistent.');
                        end
                    end
                end
                
                sizeAInCatDim = sizeA(dim);
                for ii = 1:numel(varargin)
                    subs = cell(1, ndimsA);
                    subs(:) = {':'};
                    sizeVararginInCatDim = sizeVararginElements{ii}(dim);
                    subs{dim} = sizeAInCatDim+1:sizeAInCatDim+sizeVararginInCatDim;
                    c = subsasgn(c, substruct('()', subs), varargin{ii});
                    
                    sizeAInCatDim = sizeAInCatDim+sizeVararginInCatDim;
                end
                
            end
        end
        function c = horzcat(a, varargin)
            c = cat(2, a, varargin{:});
        end
        function c = vertcat(a, varargin)
            c = cat(1, a, varargin{:});
        end
        function d = get.Value(obj)
            d = get_value(obj);
        end
        function d = get.StdUnc(obj)
            d = get_stdunc(obj);
        end
        function b = get.IsComplex(obj)
            b = LinProp.IsComplexNet(obj.NetObject);
        end
        function b = get.IsArray(obj)
            b = LinProp.IsArrayNet(obj.NetObject);
        end
        function d = double(obj)
            d = get_value(obj);
        end
        function o = get_net_object(obj)
            o = obj.NetObject;
        end
        function d = get_value(obj)
            d = LinProp.Convert2Double(LinProp.UncHelper.GetValue(obj.NetObject));
        end
        function d = get_stdunc(obj)
            d = LinProp.Convert2Double(LinProp.UncHelper.GetStdUnc(obj.NetObject));
        end
        function d = get_idof(obj)
            % GET_IDOF Inverse degree of freedom.
            %
            % d = GET_IDOF(unc) returns a matrix of the same size of unc, containing
            % the inverse degrees of freedom of every element of unc.
            d = LinProp.Convert2Double(LinProp.UncHelper.GetIDof(obj.NetObject));
        end
        function d = get_fcn_value(obj)
            d = LinProp.Convert2Double(LinProp.UncHelper.GetFcnValue(obj.NetObject));
        end
        function d = get_coverage_interval(obj, p)
            % GET_COVERAGE_INTERVAL Coverage interval bounds.
            %
            % I = get_coverage_interval(unc, p) returns a matrix of size n-by-2
            % containing the bounds of the coverage interval of unc for the probability
            % p, with 0 < p < 1 and n=numel(unc). The first column are the lower
            % bounds, the second column the upper bounds.
            %
            % The input argument unc is always interpreted as a vector, thus
            % get_coverage_interval(unc) is the same as get_coverage_interval(unc(:)).
            % If unc contains complex values, the real and imaginary parts are
            % interleaved with each imaginary part right after its respective real
            % part. This is the same as passing [real(unc(:)), imag(unc(:))]' instead
            % of unc.
            l = ToUncList(obj);
            temp = LinProp.UncHelper.GetCoverageInterval(l, p);
            array = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
            array.Init2dData(temp);
            d = LinProp.Convert2Double(array);
        end
        function d = get_moment(obj, n)
            % GET_MOMENT Central moment.
            %
            % d = GET_MOMENT(unc, n) returns the n'th central moment of the
            % distribution of obj. unc must be a scalar LinProp, n must be a
            % non-negative integer.
            d = LinProp.Convert2Double(LinProp.UncHelper.GetMoment(obj.NetObject, int32(n)));
        end
        function c = get_correlation(obj)
            % GET_CORRELATION Correlation matrix.
            %
            % C = get_correlation(unc) returns a matrix of size n-by-n containing the
            % correlation factors between the elements of unc, with n = numel(unc).
            %
            % The input argument unc is always interpreted as a vector, thus
            % get_correlation(unc) is the same as get_correlation(unc(:)). If unc
            % contains complex values, the real and imaginary parts are interleaved
            % with each imaginary part right after its respective real part. This is
            % the same as passing [real(unc(:)), imag(unc(:))]' instead of unc.
            l = ToUncList(obj);
            temp = LinProp.UncHelper.GetCorrelation(l);
            array = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
            array.Init2dData(temp);
            c = LinProp.Convert2Double(array);
        end
        function c = get_covariance(obj)
            % GET_COVARIANCE Covariance matrix.
            %
            % C = get_covariance(unc) returns a matrix of size n-by-n containing the
            % covariances of the elements of unc, with n = numel(unc).
            %
            % The input argument unc is always interpreted as a vector, thus
            % get_covariance(unc) is the same as get_covariance(unc(:)). If unc
            % contains complex values, the real and imaginary parts are interleaved
            % with each imaginary part right after its respective real part. This is
            % the same as passing [real(unc(:)), imag(unc(:))]' instead of unc.
            l = ToUncList(obj);
            temp = LinProp.UncHelper.GetCovariance(l);
            array = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
            array.Init2dData(temp);
            c = LinProp.Convert2Double(array);
        end
        function c = get_jacobi(obj)
            % GET_JACOBI Uncertainty contributions of base inputs.
            %
            % c = GET_JACOBI(unc) returns a m-by-n matrix containing the sensitivity
            % coefficients multiplied with the standard uncertainties of the base
            % inputs of unc. m is the number of elements in unc (*2 for complex
            % numbers), n is the number of base inputs of unc.
            %
            % If the order of the base inputs is not clear, use get_unc_component() and
            % explicitly specify the order of the base inputs by passing them as
            % intermediate results.
            %
            % The input argument unc is always interpreted as a vector, thus
            % get_jacobi(unc) is the same as get_jacobi(unc(:)). If unc contains
            % complex values, the real and imaginary parts are interleaved with each
            % imaginary part right after its respective real part. This is the same as
            % passing [real(unc(:)), imag(unc(:))]' instead of unc.
            % 
            % The computational and mathematical relations between base inputs and 
            % uncertainty numbers are described in <a href="https://doi.org/10.1088/0026-1394/49/6/809">doi.org/10.1088/0026-1394/49/6/809</a>.
            %
            % See also get_jacobi2, get_unc_component.
            l = ToUncList(obj);
            temp = LinProp.UncHelper.GetJacobi(l);
            array = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
            array.Init2dData(temp);
            c = LinProp.Convert2Double(array);
        end
        function c = get_jacobi2(x, y)
            % GET_JACOBI2 Sensitivities to intermediate results.
            %
            % c = GET_JACOBI2(x, y) returns a m-by-n matrix containing the sensitivities
            % of x to the intermediate results y. m is the number of elements in x (*2
            % for complex numbers), n is the number intermediate results y.
            %
            % The intermediate results y must form a base for x, thus 
            %   1. y must be complete in the sense that there exists a function f such 
            %      that x = f(y).
            %   2. the elements of y must be linearly independent, i.e. there must not
            %      exist a function g such that x = g(y') with y' being a subset of y.
            %
            % The user is required to ensure that these conditions are met.
            % Otherwise, the result returned by this method might be incorrect.
            %
            % To check that these conditions were met in the statement
            %   j = get_jacobi2(x,y)
            % test if standard uncertainty of y matches
            %   get_stdunc(x) == sqrt(j*get_covariance(y)*j')
            %
            % The input arguments x and y are always interpreted as vectors, thus
            % get_jacobi2(x, y) is the same as get_jacobi2(x(:), y(:)). If x or y
            % contain complex values, the real and imag part are treated separately,
            % thus get_jacobi2(cx, cy) is the same as get_jacobi2([real(cx(:)),
            % imag(cx(:))], [real(cy(:)), imag(cy(:))]).
            % 
            % The mathematical relations between base inputs, intermediate results, and 
            % uncertainty numbers are described in <a href="https://doi.org/10.1088/0026-1394/49/6/809">doi.org/10.1088/0026-1394/49/6/809</a>.
            %
            % See also get_jacobi, get_unc_component.
            x2 = ToUncList(x);
            y2 = ToUncList(y);
            temp = LinProp.UncHelper.GetJacobi2(x2, y2);
            array = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
            array.Init2dData(temp);
            c = LinProp.Convert2Double(array);
        end
        function c = get_unc_component(x, y)
            % GET_UNC_COMPONENT Uncertainty contributions of intermediate results.
            %
            % c = get_unc_component(x, y) returns the uncertainty contributions of the
            % intermediate results y to the uncertainty object x.
            %
            % get_unc_component(x,y) is the same as get_jacobi2(x,y) .* get_stdunc(y).
            %
            % The input arguments x and y are always interpreted as vectors, thus
            % get_unc_component(x, y) is the same as get_unc_component(x(:), y(:)). If
            % x or y contain complex values, the real and imag part are treated
            % separately, thus get_unc_component(cx, cy) is the same as
            % get_unc_component([real(cx(:)),imag(cx(:))], [real(cy(:)),imag(cy(:))]).
            % 
            % The mathematical relations between base inputs, intermediate results, and 
            % uncertainty numbers are described in <a href="https://doi.org/10.1088/0026-1394/49/6/809">doi.org/10.1088/0026-1394/49/6/809</a>.
            %
            % See also get_jacobi, get_jacobi2.
            x2 = ToUncList(x);
            y2 = ToUncList(y);
            temp = LinProp.UncHelper.GetUncComponent(x2, y2);
            array = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
            array.Init2dData(temp);
            c = LinProp.Convert2Double(array);
        end
        function n = memsize(obj)
            n = double(obj.NetObject.memsize);
        end
        function y = uplus(x)
            y = copy(x);
        end
        function y = uminus(x)
            y = LinProp(x.NetObject.Negative());
        end
        function z = plus(x,y)
            if numel(x) == 0 && numel(y) == 0
                z = LinProp([]);
            else
                x = LinProp(x);
                y = LinProp(y);
                if x.IsComplex && ~y.IsComplex
                    y = complex(y);
                end
                if ~x.IsComplex && y.IsComplex
                    x = complex(x);
                end
                if ~x.IsArray && ~y.IsArray
                    z = LinProp(x.NetObject.Add(y.NetObject));
                elseif x.IsArray && ~y.IsArray
                    z = LinProp(x.NetObject.LAdd(y.NetObject));
                elseif ~x.IsArray && y.IsArray
                    z = LinProp(y.NetObject.RAdd(x.NetObject));
                else
                    [x, y] = LinProp.replicateSingletonDimensions(x, y);
                    z = LinProp(x.NetObject.Add(y.NetObject));
                end
            end
        end
        function z = minus(x,y)
            if numel(x) == 0 && numel(y) == 0
                z = LinProp([]);
            else
                x = LinProp(x);
                y = LinProp(y);
                if x.IsComplex && ~y.IsComplex
                    y = complex(y);
                end
                if ~x.IsComplex && y.IsComplex
                    x = complex(x);
                end
                if ~x.IsArray && ~y.IsArray
                    z = LinProp(x.NetObject.Subtract(y.NetObject));
                elseif x.IsArray && ~y.IsArray
                    z = LinProp(x.NetObject.LSubtract(y.NetObject));
                elseif ~x.IsArray && y.IsArray
                    z = LinProp(y.NetObject.RSubtract(x.NetObject));
                else
                    [x, y] = LinProp.replicateSingletonDimensions(x, y);
                    z = LinProp(x.NetObject.Subtract(y.NetObject));
                end
            end
        end
        function z = times(x,y)
            x = LinProp(x);
            y = LinProp(y);
            if x.IsComplex && ~y.IsComplex
                y = complex(y);
            end
            if ~x.IsComplex && y.IsComplex
                x = complex(x);
            end
            
            if ~x.IsArray && ~y.IsArray
                z = LinProp(x.NetObject.Multiply(y.NetObject));
            elseif x.IsArray && ~y.IsArray
                z = LinProp(x.NetObject.LMultiply(y.NetObject));
            elseif ~x.IsArray && y.IsArray
                z = LinProp(y.NetObject.RMultiply(x.NetObject));
            else
                [x, y] = LinProp.replicateSingletonDimensions(x, y);
                z = LinProp(x.NetObject.Multiply(y.NetObject));
            end
        end
        function z = rdivide(x,y)
            x = LinProp(x);
            y = LinProp(y);
            if x.IsComplex && ~y.IsComplex
                y = complex(y);
            end
            if ~x.IsComplex && y.IsComplex
                x = complex(x);
            end
            if ~x.IsArray && ~y.IsArray
                z = LinProp(x.NetObject.Divide(y.NetObject));
            elseif x.IsArray && ~y.IsArray
                z = LinProp(x.NetObject.LDivide(y.NetObject));
            elseif ~x.IsArray && y.IsArray
                z = LinProp(y.NetObject.RDivide(x.NetObject));
            else
                [x, y] = LinProp.replicateSingletonDimensions(x, y);
                z = LinProp(x.NetObject.Divide(y.NetObject));
            end
        end
        function z = power(x,y)
            x = LinProp(x);
            y = LinProp(y);
            ydbl = double(y);
            yint = int32(ydbl);
            if (~y.IsArray) && (~y.IsComplex)
                yconst = (yint == ydbl) & (y.NetObject.IsConst);
            else
                yconst = false;
            end
            if yconst
                if x.IsArray
                    z = LinProp.Convert2LinProp(x.NetObject.Pow(yint));
                else
                    z = LinProp.Convert2LinProp(x.NetObject.Pow(yint));
                end
            else
                value = get_value(x); 
                if any(value(:) < 0)
                    x = complex(x);
                end
                if x.IsComplex && ~y.IsComplex
                    y = complex(y);
                end
                if ~x.IsComplex && y.IsComplex
                    x = complex(x);
                end
                if ~x.IsArray && ~y.IsArray
                    z = LinProp.Convert2LinProp(x.NetObject.Pow(y.NetObject));
                elseif x.IsArray && ~y.IsArray
                    y = y.*ones(size(x));
                    z = x.^y;
                    % z = LinProp.Convert2LinProp(x.NetObject.LPow(y.NetObject));
                elseif ~x.IsArray && y.IsArray
                    x = x.*ones(size(y));
                    z = x.^y;
                    % z = LinProp.Convert2LinProp(y.NetObject.RPow(x.NetObject));
                else
                    z = LinProp.Convert2LinProp(x.NetObject.Pow(y.NetObject));
                end
            end
        end
        function z = complex(varargin)
            narginchk(1, 2);
            if nargin == 1
                x = varargin{1};
                if x.IsComplex
                    z = copy(x);
                else
                    if x.IsArray
                        z = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    else
                        z = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                    end
                    z.InitRe(x.NetObject);
                    z = LinProp(z);
                end
            else
                a = LinProp(varargin{1});
                b = LinProp(varargin{2});
                if ~isreal(a)
                    error('Input for real part must be a real-valued.');
                end
                if ~isreal(b)
                    error('Input for imaginary part must be a real-valued.');
                end
                if ~isscalar(a) && ~isscalar(b)
                    if ~isequal(size(a), size(b))
                        throwAsCaller(MException('MATLAB:sizeDimensionsMustMatch', 'Input arrays must have the same size.'));
                    end
                else
                    [a, b] = LinProp.replicateSingletonDimensions(a, b);
                end
                if a.IsArray
                    z = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                else
                    z = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                end
                z.InitReIm(a.NetObject, b.NetObject);
                z = LinProp(z);
            end
        end
        function y = real(x)
            if x.IsComplex
                y = LinProp(x.NetObject.Real());
            else
                y = copy(x);
            end
        end
        function y = imag(x)
            x = complex(x);
            y = LinProp(x.NetObject.Imag());
        end
        function y = conj(x)
            x = complex(x);
            y = LinProp(x.NetObject.Conj());
        end        
        function y = abs(x)
            y = LinProp(x.NetObject.Abs());
        end       
        function y = angle(x)
            x = complex(x);
            y = LinProp(x.NetObject.Angle());
        end
        function q = unwrap(p, varargin)
            q = p + unwrap(double(p), varargin{:}) - double(p);
        end
        function y = deg2rad(x)
            y = (pi/180) .* x;
        end
        function y = rad2deg(x)
            y = (180/pi) .* x;
        end
        function y = exp(x)
            y = LinProp(x.NetObject.Exp());
        end
        function y = log(x)
            value = get_value(x); 
            if any(value(:) < 0)
                x = complex(x);
            end
            y = LinProp(x.NetObject.Log());
        end
        function y = log10(x)
            value = get_value(x); 
            if any(value(:) < 0)
                x = complex(x);
            end
            y = LinProp(x.NetObject.Log10());
        end
        function y = sqrt(x)
            value = get_value(x); 
            if any(value(:) < 0)
                x = complex(x);
            end
            y = LinProp(x.NetObject.Sqrt());
        end
        function y = sign(x)
            y = sign(double(x));
        end
        function y = sin(x)
            y = LinProp(x.NetObject.Sin());
        end
        function y = cos(x)
            y = LinProp(x.NetObject.Cos());
        end
        function y = tan(x)
            y = LinProp(x.NetObject.Tan());
        end
        function y = sinh(x)
            y = LinProp(x.NetObject.Sinh());
        end
        function y = cosh(x)
            y = LinProp(x.NetObject.Cosh());
        end
        function y = tanh(x)
            y = LinProp(x.NetObject.Tanh());
        end
        function y = asin(x)
            y = LinProp(x.NetObject.Asin());
        end
        function y = acos(x)
            y = LinProp(x.NetObject.Acos());
        end
        function y = atan(x)
            y = LinProp(x.NetObject.Atan());
        end
        function z = atan2(x,y)
            x = LinProp(x);
            y = LinProp(y);
            if x.IsComplex || y.IsComplex
                error('Inputs must be real');
            end
            if ~x.IsArray && ~y.IsArray
                z = LinProp(x.NetObject.Atan2(y.NetObject));
            elseif x.IsArray && ~y.IsArray
                y = y.*ones(size(x));
                z = atan2(x,y);
                % z = LinProp(x.NetObject.LAtan2(y.NetObject));
            elseif ~x.IsArray && y.IsArray
                x = x.*ones(size(y));
                z = atan2(x,y);
                % z = LinProp(y.NetObject.RAtan2(x.NetObject));
            else
                z = LinProp(x.NetObject.Atan2(y.NetObject));
            end
        end
        function y = asinh(x)
            x = complex(x);
            y = LinProp(x.NetObject.Asinh());
        end
        function y = acosh(x)
            x = complex(x);
            y = LinProp(x.NetObject.Acosh());
        end
        function y = atanh(x)
            x = complex(x);
            y = LinProp(x.NetObject.Atanh());
        end
        function [k,e] = ellipke(x)
            if (x.IsComplex)
                error('Input must be real');
            end
            k = LinProp(x.NetObject.Ellipk());
            e = LinProp(x.NetObject.Ellipe());
        end
        function z = eq(x,y)
            z = double(x) == double(y);
        end
        function z = ge(x,y)
            z = double(x) >= double(y);
        end
        function z = gt(x,y)
            z = double(x) > double(y);
        end
        function z = le(x,y)
            z = double(x) <= double(y);
        end
        function z = lt(x,y)
            z = double(x) < double(y);
        end
        function z = ne(x,y)
            z = double(x) ~= double(y);
        end
        function y = isfinite(x)
            y = isfinite(double(x));
        end
        function y = isinf(x)
            y = isinf(double(x));
        end
        function y = isnan(x)
            y = isnan(double(x));
        end
        function y = isreal(x)
            y = ~x.IsComplex;
        end
        function y = transpose(x)
            if x.IsArray
                y = LinProp(x.NetObject.Transpose());
            else
                y = x;
            end
        end
        function y = ctranspose(x)
            if x.IsArray
                if x.IsComplex
                    y = LinProp(x.NetObject.CTranspose());
                else
                    y = LinProp(x.NetObject.Transpose());
                end
            else
                y = conj(x);
            end
        end
        function d = diag(A)
            am = LinProp.Convert2UncArray(A);
            s = size(A);
            if ((s(1) == 1) || (s(2) == 1))
                n = numel(A);
                if A.IsComplex
                    m = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    m.Zeros2d(n, n);
                else
                    m = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    m.Zeros2d(n, n);
                end
                for i1 = 1:n
                    m.SetItem2d(i1-1, i1-1, am.GetItem1d(i1-1));
                end
                d = LinProp.Convert2LinProp(m);
            else
                if ((am.ndims ~= 2) || (s(1) ~= s(2)))
                    error('Matrix must be square.');
                end
                n1 = s(1);
                if A.IsComplex
                    m = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    m.Init2d(n1, 1);
                else
                    m = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    m.Init2d(n1, 1);
                end
                for i1 = 1:n1
                    m.SetItem2d(i1-1, 0, am.GetItem2d(i1-1, i1-1));
                end
                d = LinProp.Convert2LinProp(m);
            end
        end
        function d = det(A)
            linalg = LinProp.LinAlg(A.IsComplex);
            d = LinProp.Convert2LinProp(linalg.Det(LinProp.Convert2UncArray(A)));
        end
        function d = inv(A)
            linalg = LinProp.LinAlg(A.IsComplex);
            d = LinProp.Convert2LinProp(linalg.Inv(LinProp.Convert2UncArray(A)));
        end
        function z = ldivide(x,y)
            z = y./x;
        end
        function z = mldivide(x,y)
            x = LinProp(x);
            y = LinProp(y);
            if size(y, 2) == 1
                if x.IsComplex && ~y.IsComplex
                    y = complex(y);
                end
                if ~x.IsComplex && y.IsComplex
                    x = complex(x);
                end
                xm = LinProp.Convert2UncArray(x);
                yv = LinProp.Convert2UncArray(y);
                s = size(x);
                if s(1) == s(2)
                    linalg = LinProp.LinAlg(x.IsComplex);
                    zv = linalg.Solve(xm, yv);
                else
                    linalg = LinProp.LinAlg2(x.IsComplex);
                    zv = linalg.LstSqrSolve(xm, yv);
                end
                z = LinProp.Convert2LinProp(zv);
            else
                z = inv(x)*y;
            end
        end
        function z = mpower(x,y)
            if (numel(x) == 1 && numel(y) == 1)
                % Power
                z = x.^y;
            elseif (numel(x) == 1)
                % Matrix Exponents
                [v, d] = eig(y);
                z = v*diag(x.^diag(d))/(v);
            elseif (numel(y) == 1)
                if (isa(y, 'double') && y == -1)
                    % Matrix Inverse
                    z = inv(x);
                else
                    % Matrix Power
                    [v, d] = eig(x);
                    z = v*diag(diag(d).^y)/(v);
                end
            else
                error('Incorrect dimensions for raising a matrix to a power. Check that the matrix is square and the power is a scalar. To perform elementwise matrix powers, use ''.^''.')
            end
        end
        function z = mrdivide(x,y)
            z = x*inv(y);
        end
        function z = mtimes(x,y)
            if isscalar(x) || isscalar(y)
                z = times(x, y);
            else
                x = LinProp(x);
                y = LinProp(y);
                if x.IsComplex && ~y.IsComplex
                    y = complex(y);
                end
                if ~x.IsComplex && y.IsComplex
                    x = complex(x);
                end
                
                dims = max(ndims(x), ndims(y));
                if dims > 2
                    error('Arguments must be 2-D, or at least one argument must be scalar. Use TIMES (.*) for elementwise multiplication.');
                elseif size(x, 2) ~= size(y, 1)
                    error('Incorrect dimensions for matrix multiplication. Check that the number of columns in the first matrix matches the number of rows in the second matrix. To perform elementwise multiplication, use ''.*''.');
                end
                
                linalg = LinProp.LinAlg(x.IsComplex);
                xm = LinProp.Convert2UncArray(x);
                ym = LinProp.Convert2UncArray(y);
                zm = linalg.Dot(xm, ym);
                z = LinProp.Convert2LinProp(zm);
            end
        end
        function [L, U, P] = lu(A)
            linalg = LinProp.LinAlg(A.IsComplex);
            am = LinProp.Convert2UncArray(A);
            temp = linalg.Lu(am);
            L = LinProp.Convert2LinProp(temp.l);
            U = LinProp.Convert2LinProp(temp.u);
            P = LinProp.Convert2LinProp(temp.p);
        end
        function U = chol(A)
            linalg = LinProp.LinAlg2(A.IsComplex);
            am = LinProp.Convert2UncArray(A);
            lm = linalg.Cholesky(am);
            L = LinProp.Convert2LinProp(lm);
            U = L';
        end
        function [Q, R] = qr(A)
            linalg = LinProp.LinAlg2(A.IsComplex);
            am = LinProp.Convert2UncArray(A);
            [qm, rm] = linalg.Qr(am);
            Q = LinProp.Convert2LinProp(qm);
            R = LinProp.Convert2LinProp(rm);
        end
        function [U, S, V] = svd(A)
            linalg = LinProp.LinAlg2(A.IsComplex);
            am = LinProp.Convert2UncArray(A);
            [um, sm, vm] = linalg.Svd(am);
            U = LinProp.Convert2LinProp(um);
            S = LinProp.Convert2LinProp(sm);
            V = LinProp.Convert2LinProp(vm);
        end
        function x = lscov(A,b,V)
            A = LinProp(A);
            b = LinProp(b);
            V = LinProp(V);
            if A.IsComplex && ~b.IsComplex
                b = complex(b);
            end
            if ~A.IsComplex && b.IsComplex
                A = complex(A);
            end
            linalg = LinProp.LinAlg2(A.IsComplex);
            Am = LinProp.Convert2UncArray(A);
            bv = LinProp.Convert2UncArray(b);
            if isvector(V)
                W = diag(V);
                Wm = LinProp.Convert2UncArray(W);
                xv = linalg.WeightedLstSqrSolve(Am, bv, Wm);
            else
                Vm = LinProp.Convert2UncArray(V);
                xv = linalg.GeneralLstSqrSolve(Am, bv, Vm);
            end
            x = LinProp.Convert2LinProp(xv);
        end
        function a = sum(x, varargin)
            s = size(x);
            switch nargin
                case 1
                    % find first non-singleton dimension
                    f = [find(s > 1) 1];
                    a = sum(x, f(1));
                case {2,3}
                    i = varargin{1};
                    if i < 1
                        error('Dimension argument must be a positive integer scalar within indexing range');
                    end
                    if i > numel(s)
                        a = x;
                    else
                        linalg = LinProp.LinAlg(x.IsComplex);
                        xm = LinProp.Convert2UncArray(x);
                        am = linalg.Sum(xm, i-1);
                        a = LinProp.Convert2LinProp(am);
                    end
                otherwise
                    error('Too many input arguments')
            end
        end
        function a = prod(x, varargin)
            s = size(x);
            switch nargin
                case 1
                    % find first non-singleton dimension
                    f = [find(s > 1) 1];
                    a = prod(x, f(1));
                case {2,3}
                    i = varargin{1};
                    if i < 1
                        error('Dimension argument must be a positive integer scalar within indexing range');
                    end
                    if i > numel(s)
                        a = x;
                    else
                        linalg = LinProp.LinAlg(x.IsComplex);
                        xm = LinProp.Convert2UncArray(x);
                        am = linalg.Prod(xm, i-1);
                        a = LinProp.Convert2LinProp(am);
                    end
                otherwise
                    error('Too many input arguments')
            end
        end
        function a = mean(x, varargin)
            s = size(x);
            switch nargin
                case 1
                    % find first non-singleton dimension
                    f = [find(s > 1) 1];
                    a = mean(x, f(1));
                case 2
                    i = varargin{1};
                    if i < 1
                        error('Dimension argument must be a positive integer scalar within indexing range');
                    end
                    if i > numel(s)
                        a = x;
                    else
                        a = sum(x, i)./size(x, i);
                    end
                otherwise
                    error('Too many input arguments')
            end
        end
        function X = fft(A)
            numlib = LinProp.NumLib2(1);
            A = complex(A);
            s = size(A);
            am = LinProp.Convert2UncArray(A);
            xm = numlib.Fft(am);
            X = LinProp.Convert2LinProp(xm);
            X = reshape(X, s);
        end
        function X = ifft(A)
            numlib = LinProp.NumLib2(1);
            A = complex(A);
            s = size(A);
            am = LinProp.Convert2UncArray(A);
            xm = numlib.Ifft(am);
            X = LinProp.Convert2LinProp(xm);
            X = reshape(X, s);
        end
        function yy = interpolation(x, y, n, xx)
            % yy = INTERPOLATION(x, y, n, xx) Interpolation.
            %
            % Interpolates y(x) at points xx, with y(x) being a polynomial of n-th
            % degree, specified by the vectors x and y. Returns yy, a LinProp vector of
            % the same size as xx, which contains the interpolated values. The
            % uncertainties of y (i.e. y(x)) are propagated, while any uncertainties of
            % x and xx are ignored. While y must be a LinProp, x and xx can be any
            % type. The parameter n must be a positive integer smaller than numel(x).
            % 
            % In general, the interpolated values are calculated based on n values of x
            % and y. Using interpolation, this will result in the uncertainties of yy
            % being smaller (or at the edges larger) that those of y. Using interpolation2, 
            % the uncertainties of yy will be a linear interpolation of y.
            % See <a href="matlab:s=which('LinProp');[s,~,~]=fileparts(s);edit([s,'\..\Examples\Example_Interpolation.m']);">Examples/Example_Interpolation.m</a>
            %
            % See also LinProp.interpolation2, LinProp.spline.
            x = double(x(:));
            y = LinProp(y);
            n = int32(n);
            s = size(xx);
            xx = double(xx(:));
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            yym = numlib.Interpolation(x, ym, n, xx);
            yy = LinProp.Convert2LinProp(yym);
            yy = reshape(yy, s);
        end
        function yy = interpolation2(x, y, n, xx)
            % yy = INTERPOLATION2(x, y, n, xx) Interpolation with linear unc. propagation.
            %
            % Interpolates y(x) at points xx, with y(x) being a polynomial of n-th
            % degree, specified by the vectors x and y. Returns yy, a LinProp vector of
            % the same size as xx, which contains the interpolated values. The
            % uncertainties of y (i.e. y(x)) are linearly interpolated, while any
            % uncertainties of x and xx are ignored. While y must be a LinProp, x and
            % xx can be any type. The parameter n must be a positive integer smaller
            % than numel(x).
            % 
            % In general, the interpolated values are calculated based on n values of x
            % and y. Using interpolation, this will result in the uncertainties of yy
            % being smaller (or at the edges larger) that those of y. Using interpolation2, 
            % the uncertainties of yy will be a linear interpolation of y.
            % See <a href="matlab:s=which('LinProp');[s,~,~]=fileparts(s);edit([s,'\..\Examples\Example_Interpolation.m']);">Examples/Example_Interpolation.m</a>
            %
            % See also LinProp.interpolation, LinProp.spline2.
            x = double(x(:));
            y = LinProp(y);
            n = int32(n);
            s = size(xx);
            xx = double(xx(:));
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            yym = numlib.Interpolation2(x, ym, n, xx);
            yy = LinProp.Convert2LinProp(yym);
            yy = reshape(yy, s);
        end
        function yy = spline(x, y, xx, varargin)
            % yy = SPLINE(x, y, xx, [bounds]) Spline interpolation.
            %
            % Interpolates y(x) at points xx, with y(x) being a cubic spline defined by
            % the vectors x and y and the boundary conditions. Returns yy, a LinProp
            % vector of the same size as xx, which contains the interpolated values.
            % The uncertainties of y(x) are propagated, while any uncertainties of x
            % and xx are ignored. While y must be a LinProp, x and xx can be any type.
            %
            % SPLINE(x, y) Uses the 'not-a-knot' boundary condition.
            %
            % SPLINE(__, boundaryCond) Specifies the boundary condition on both ends.
            % Valid values for boundaryCond are 'not-a-knot' (default),
            % 'natural spline', '1st derivative', and '2nd derivative'. With the 
            % conditions '1st derivative' and '2nd derivative', the respective 
            % derivatives are set to zero.
            % 
            % SPLINE(__, leftBoundCond, leftValue, rightBoundCond, rightValue)
            % Specifies the boundary condition for the left and right end and also the
            % value of the derivatives. Valid values for leftBoundCond and
            % rightBoundCond are 'not-a-knot' (default), 'natural spline', 
            % '1st derivative', and '2nd derivative'.
            %
            % In general, the interpolated values are calculated based on multiple
            % values of x and y. Using spline, this will result in the uncertainties of
            % yy being smaller (or at the edges larger) that those of y. Using spline2,
            % the uncertainties of yy will be a linear interpolation of y. See 
            % <a href="matlab:s=which('LinProp');[s,~,~]=fileparts(s);edit([s,'\..\Examples\Example_Interpolation.m']);">Examples/Example_Interpolation.m</a>
            %
            % See also LinProp.interpolation, LinProp.spline2, LinProp.splinecoefs.
            x = double(x(:));
            y = LinProp(y);
            s = size(xx);
            xx = double(xx(:));
            [y, sb, sv, eb, ev] = SplineOptArgs(y, varargin{:});
            numlib = LinProp.NumLib(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            yym = numlib.SplineInterpolation(x, ym, xx, sb, sv, eb, ev);
            yy = LinProp.Convert2LinProp(yym);
            yy = reshape(yy, s);
        end
        function yy = spline2(x, y, xx, varargin)
            % yy = SPLINE2(x, y, xx, [bounds]) Spline interpolation with linear unc. propagation.
            %
            % Interpolates y(x) at points xx, with y(x) being a cubic spline defined by
            % the vectors x and y and the boundary conditions. Returns yy, a LinProp
            % vector of the same size as xx, which contains the interpolated values. The
            % uncertainties of y(x) are linearly propagated, while any uncertainties of
            % x and xx are ignored. While y must be a LinProp, x and xx can be any
            % type.
            %
            % SPLINE2(x, y) Uses the 'not-a-knot' boundary condition.
            %
            % SPLINE2(__, boundaryCond) Specifies the boundary condition on both ends.
            % Valid values for boundaryCond are 'not-a-knot' (default),
            % 'natural spline', '1st derivative', and '2nd derivative'. With the 
            % conditions '1st derivative' and '2nd derivative', the respective 
            % derivatives are set to zero.
            % 
            % SPLINE2(__, leftBoundCond, leftValue, rightBoundCond, rightValue)
            % Specifies the boundary condition for the left and right end and also the
            % value of the derivatives. Valid values for leftBoundCond and
            % rightBoundCond are 'not-a-knot' (default), 'natural spline', 
            % '1st derivative', and '2nd derivative'.
            %
            % In general, the interpolated values are calculated based on multiple
            % values of x and y. Using spline, this will result in the uncertainties of
            % yy being smaller (or at the edges larger) that those of y. Using spline2,
            % the uncertainties of yy will be a linear interpolation of y. See 
            % <a href="matlab:s=which('LinProp');[s,~,~]=fileparts(s);edit([s,'\..\Examples\Example_Interpolation.m']);">Examples/Example_Interpolation.m</a>
            %
            % See also LinProp.interpolation2, LinProp.spline, LinProp.splinecoefs.
            x = double(x(:));
            y = LinProp(y);
            s = size(xx);
            xx = double(xx(:));
            [y, sb, sv, eb, ev] = SplineOptArgs(y, varargin{:});
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            yym = numlib.SplineInterpolation2(x, ym, xx, sb, sv, eb, ev);
            yy = LinProp.Convert2LinProp(yym);
            yy = reshape(yy, s);
        end
        function p = splinecoefs(x, y, varargin)
            % p = SPLINECOEFS(x, y, [bounds]) Coefficients of interpolation spline.
            %
            % Returns the coefficients of the cubic splines that connect the points
            % defined by the vectors x and y and the specified the boundary conditions.
            % Returns a n-by-4 matrix of local coefficients, where n is the number of
            % spline segments, i.e. n = length(x)-1. While x and y can be any type,
            % uncertainties of x are ignored.
            %
            % The i'th row of p with values [a b c d] can be interpreted as 
            %   f(xq) = a*(xq - x(i)).^3 + b*(xq - x(i)).^2 + c*(xq - x(i)) + d.
            % The matrix returned by this method has the same properties as the
            % pp.coefs matrix returned by <a href="matlab:help spline">PP = spline(X,Y)</a>.
            %
            % SPLINECOEFS(x, y) Uses the 'not-a-knot' boundary condition.
            %
            % SPLINECOEFS(__, boundaryCond) Specifies the boundary condition on both
            % ends. Valid values for boundaryCond are 'not-a-knot' (default),
            % 'natural spline', '1st derivative', and '2nd derivative'. With the 
            % conditions '1st derivative' and '2nd derivative', the respective 
            % derivatives are set to zero.
            % 
            % SPLINECOEFS(__, leftBoundCond, leftValue, rightBoundCond, rightValue)
            % Specifies the boundary condition for the left and right end and also the
            % value of the derivatives. Valid values for leftBoundCond and
            % rightBoundCond are 'not-a-knot' (default), 'natural spline', 
            % '1st derivative', and '2nd derivative'.
            %
            % See also LinProp.spline.
            x = double(x(:));
            y = LinProp(y);
            [y, sb, sv, eb, ev] = SplineOptArgs(y, varargin{:});
            numlib = LinProp.NumLib(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            pm = numlib.SplineCoefs(x, ym, sb, sv, eb, ev);
            p = LinProp.Convert2LinProp(pm);
        end
        function a = integrate(x, y, n)
            % a = INTEGRATE(x, y, n) Integration with cumulative result.
            %
            % Calculates the numerical integral of y(x), with y(x) being a polynomial
            % of n-th degree specified by the vectors x and y. Returns a LinProp vector
            % of the same size as y which contains the cumulative integral up to every
            % value of x. The uncertainties of y(x) are propagated, while any uncer-
            % tainties of x are ignored. While y must be a LinProp, x can be any
            % type. The parameter n must be a positive integer smaller than numel(y).
            %
            % See also LinProp.integrate2, LinProp.splineintegrate.

            x = double(x(:));
            y = LinProp(y);
            n = int32(n);
            s = size(y);
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            am = numlib.Integrate(x, ym, n);
            a = LinProp.Convert2LinProp(am);
            a = reshape(a, s);
        end
        function a = integrate2(x, y, n)
            % a = INTEGRATE2(x, y, n) Integration with scalar result.
            %
            % Calculates the numerical integral of y(x), with y(x) being a polynomial
            % of n-th degree specified by the vectors x and y. Returns the result of
            % the whole integral as a LinProp scalar. The input arguments x and y
            % specify y(x). The uncertainties of y(x) are propagated, while any
            % uncertainties of x are ignored. While y must be a LinProp, x can be any
            % type. The parameter n must be a positive integer smaller than numel(y).
            %
            % a = integrate2(x, y, n) returns the same result as:
            %   a = integrate(x, y, n);
            %   a = a(end);
            %
            % See also LinProp.integrate, LinProp.splineintegrate2.
            
            x = double(x(:));
            y = LinProp(y);
            n = int32(n);
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            am = numlib.Integrate2(x, ym, n);
            a = LinProp.Convert2LinProp(am);
        end
        function a = splineintegrate(x, y, varargin)
            % a = SPLINEINTEGRATE(x, y, ...) Spline integration with cumulative result.
            %
            % Calculates the numerical integral of y(x), with y(x) being a cubic spline
            % defined by the vectors x and y and the boundary conditions. Returns a
            % LinProp vector of the same size as y which contains the cumulative
            % integral up to every value of x. The uncertainties of y(x) are
            % propagated, while any uncertainties of x are ignored. While y has to be a
            % LinProp, x can be any type.
            %
            % SPLINEINTEGRATE(x, y) Uses the 'not-a-knot' boundary condition.
            %
            % SPLINEINTEGRATE(__, boundaryCond) Specifies the boundary condition on
            % both ends. Valid values for boundaryCond are 'not-a-knot' (default),
            % 'natural spline', '1st derivative', and '2nd derivative'. With the
            % conditions '1st derivative' and '2nd derivative', the respective
            % derivatives are set to zero.
            % 
            % SPLINEINTEGRATE(__, leftBoundCond, leftValue, rightBoundCond, rightValue)
            % Specifies the boundary condition for the left and right end and also the
            % value of the derivatives. Valid values for leftBoundCond and
            % rightBoundCond are 'not-a-knot' (default), 'natural spline', 
            % '1st derivative', and '2nd derivative'.
            %
            % See also LinProp.integrate, LinProp.splineintegrate2.
            x = double(x(:));
            y = LinProp(y);
            s = size(y);
            [y, sb, sv, eb, ev] = SplineOptArgs(y, varargin{:});
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            am = numlib.SplineIntegrate(x, ym, sb, sv, eb, ev);
            a = LinProp.Convert2LinProp(am);
            a = reshape(a, s);
        end
        function a = splineintegrate2(x, y, varargin)
            % a = SPLINEINTEGRATE2(x, y, ...) Spline integration with scalar result.
            %
            % Calculates the numerical integral of y(x), with y(x) being a cubic spline
            % defined by the vectors x and y and the boundary conditions. Returns the
            % result of the whole integral as a LinProp scalar. The uncertainties of
            % y(x) are propagated, while any uncertainties of x are ignored. While y
            % must be a LinProp, x can be any type.
            %
            % SPLINEINTEGRATE2(x, y) Uses the 'not-a-knot' boundary condition.
            %
            % SPLINEINTEGRATE2(__, boundaryCond) Specifies the boundary condition on
            % both ends. Valid values for boundaryCond are 'not-a-knot' (default),
            % 'natural spline', '1st derivative', and '2nd derivative'. With the
            % conditions '1st derivative' and '2nd derivative', the respective
            % derivatives are set to zero.
            % 
            % SPLINEINTEGRATE2(__, leftBoundCond, leftValue, rightBoundCond, rightValue)
            % Specifies the boundary condition for the left and right end and also the
            % value of the derivatives. Valid values for leftBoundCond and
            % rightBoundCond are 'not-a-knot' (default), 'natural spline', 
            % '1st derivative', and '2nd derivative'.
            %
            % See also LinProp.integrate2, LinProp.splineintegrate.
            x = double(x(:));
            y = LinProp(y);
            [y, sb, sv, eb, ev] = SplineOptArgs(y, varargin{:});
            numlib = LinProp.NumLib2(y.IsComplex);
            ym = LinProp.Convert2UncArray(y);
            am = numlib.SplineIntegrate2(x, ym, sb, sv, eb, ev);
            a = LinProp.Convert2LinProp(am);
         end
        function p = polyfit(x,y,n)
            x = LinProp(x);
            y = LinProp(y);
            n = int32(n);
            if x.IsComplex && ~y.IsComplex
                y = complex(y);
            end
            if ~x.IsComplex && y.IsComplex
                x = complex(x);
            end
            numlib = LinProp.NumLib2(x.IsComplex);
            xm = LinProp.Convert2UncArray(x);
            ym = LinProp.Convert2UncArray(y);
            pm = numlib.PolyFit(xm, ym, n);
            p = LinProp.Convert2LinProp(pm);
        end
        function y = polyval(p,x)
            p = LinProp(p);
            x = LinProp(x);
            if p.IsComplex && ~x.IsComplex
                x = complex(x);
            end
            if ~p.IsComplex && x.IsComplex
                p = complex(p);
            end
            numlib = LinProp.NumLib(p.IsComplex);
            pm = LinProp.Convert2UncArray(p);
            xm = LinProp.Convert2UncArray(x);
            ym = numlib.PolyVal(pm, xm);
            y = LinProp.Convert2LinProp(ym);
        end
        function binary_file(x, filepath)
            x.NetObject.BinarySerialize(filepath);
        end
        function xml_file(x, filepath)
            x.NetObject.XmlSerialize(filepath);
        end
        function s = xml_string(x)
            s = char(x.NetObject.XmlSerializeToString());
        end
        function bin = saveobj(obj)
            bin.data = uint8(obj.NetObject.BinarySerializeToByteArray());
            bin.array = obj.IsArray;
            bin.complex = obj.IsComplex;
        end
    end
    methods(Access = private)
        function l = ToUncList(obj)
            temp = Metas.UncLib.LinProp.UncList();
            l = temp.op_Implicit(obj.NetObject);
        end
        function [y, sb, sv, eb, ev] = SplineOptArgs(y, varargin)
            % (..., y)
            % (..., boundaryConditions)
            % (..., startBoundaryCondition, startV, endBoundaryCondition, endV)
            %
            % Boundary conditions are strings. Valid values are:
            %   'not-a-knot' (default)
            %   'natural spline'
            %   'first derivative'
            %   'second derivative'
            switch nargin
                case 1
                    sb = Metas.UncLib.Core.SplineBoundary.Not_a_Knot;
                    sv = LinProp(0);
                    eb = Metas.UncLib.Core.SplineBoundary.Not_a_Knot;
                    ev = LinProp(0);
                case 2
                    sb = BoundaryArg(varargin{1});
                    sv = LinProp(0);
                    eb = sb;
                    ev = LinProp(0);
                 case 5
                    sb = BoundaryArg(varargin{1});
                    sv = LinProp(varargin{2});
                    sv = sv(1);
                    eb = BoundaryArg(varargin{3});
                    ev = LinProp(varargin{4});
                    ev = ev(1);
                otherwise
                    error('Wrong number of input arguments')
            end
            if (y.IsComplex || sv.IsComplex || ev.IsComplex)
                y = complex(y);
                sv = complex(sv);
                ev = complex(ev);
            end
            sv = sv.NetObject;
            ev = ev.NetObject;
            
            function c = BoundaryArg(b)
                if (ischar(b))
                    b = lower(b);
                    switch b
                        case 'not-a-knot'
                            c = Metas.UncLib.Core.SplineBoundary.Not_a_Knot;
                        case 'not a knot'
                            c = Metas.UncLib.Core.SplineBoundary.Not_a_Knot;
                        case 'natural spline'
                            c = Metas.UncLib.Core.SplineBoundary.Natural_Spline;
                        case 'first derivative'
                            c = Metas.UncLib.Core.SplineBoundary.First_Derivative;
                        case '1st derivative'
                            c = Metas.UncLib.Core.SplineBoundary.First_Derivative;
                        case 'second derivative'
                            c = Metas.UncLib.Core.SplineBoundary.Second_Derivative;
                        case '2nd derivative'
                            c = Metas.UncLib.Core.SplineBoundary.Second_Derivative;
                        otherwise
                            error('Unknown spline boundary mode')
                    end
                elseif (isa(b,'Metas.UncLib.Core.SplineBoundary'))
                    c = b;
                else
                    error('Unknown spline boundary type')
                end
            end
        end
    end
    methods(Static = true)
        function obj = loadobj(bin)
            UncPropLoadNETAssemblies('LinProp');
            if bin.array
                if bin.complex
                    t = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                else
                    t = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                end
            else
                if bin.complex
                    t = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                else
                    t = Metas.UncLib.LinProp.UncNumber();
                end
            end
            v = t.BinaryDeserializeFromByteArray(bin.data(:));
            obj = LinProp(v);
        end
        % Support for array creation functions.
        % See: https://www.mathworks.com/help/releases/R2021a/matlab/matlab_oop/class-support-for-array-creation-functions.html
        function x = zeros(varargin)
            x = LinProp(zeros(varargin{:}));
        end
        function x = ones(varargin)
            x = LinProp(ones(varargin{:}));
        end
        function x = eye(varargin)
            x = LinProp(eye(varargin{:}));
        end
        function x = nan(varargin)
            x = LinProp(nan(varargin{:}));
        end
        function x = inf(varargin)
            x = LinProp(inf(varargin{:}));
        end
        function x = rand(varargin)
            x = LinProp(rand(varargin{:}));
        end
        function x = randi(varargin)
            x = LinProp(randi(varargin{:}));
        end
        function x = randn(varargin)
            x = LinProp(randn(varargin{:}));
        end
        function x = empty(varargin)
            try
                x = reshape(LinProp([]), varargin{:});
            catch e
                switch (e.identifier)
                    case 'MATLAB:getReshapeDims:notSameNumel'
                        error('MATLAB:class:emptyMustBeZero', 'At least one dimension must be zero.');
                    case 'MATLAB:getReshapeDims:sizeVector'
                        % This error is triggered with LinProp.empty(0),
                        % which should return a 0-by-0 element.
                        x = LinProp([]);
                    otherwise
                        rethrow(e);
                end
            end
        end
    end
    properties (Constant, Access = private)
        UncHelper = LinProp.UncHelperFactory();
    end
    methods(Static = true, Access = private)
        function [x, y] = replicateSingletonDimensions(x, y)
            dims = max(ndims(x), ndims(y));
            sizeX = size(x, 1:dims);
            sizeY = size(y, 1:dims);
            if any(sizeX ~= sizeY & sizeX ~= 1 & sizeY ~= 1)
                throwAsCaller(MException('MATLAB:sizeDimensionsMustMatch', 'Arrays have incompatible sizes for this operation.'));
            end
            doRepX = sizeX ~= sizeY & sizeX == 1;
            if any(doRepX)
                repX = ones(1, dims);
                repX(doRepX) = sizeY(doRepX);
                x = repmat(x, repX);
            end

            doRepY = sizeY ~= sizeX & sizeY == 1;
            if any(doRepY)
                repY = ones(1, dims);
                repY(doRepY) = sizeX(doRepY);
                y = repmat(y, repY);
            end
        end
        function h = UncHelperFactory()
            UncPropLoadNETAssemblies('LinProp');
            h = NET.createGeneric('Metas.UncLib.Core.Unc.GenericUnc', {'Metas.UncLib.LinProp.UncList', 'Metas.UncLib.LinProp.UncNumber'});
        end
        function l = LinAlg(complex)
            if complex
                lu_res = NET.GenericClass('Metas.UncLib.Core.Ndims.ComplexLuResult', 'Metas.UncLib.LinProp.UncNumber');
                narray = NET.GenericClass('Metas.UncLib.Core.Ndims.ComplexNArray', 'Metas.UncLib.LinProp.UncNumber');
                number = NET.GenericClass('Metas.UncLib.Core.Complex', 'Metas.UncLib.LinProp.UncNumber');
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.LinAlg', {lu_res, narray, number});
            else
                lu_res = NET.GenericClass('Metas.UncLib.Core.Ndims.RealLuResult', 'Metas.UncLib.LinProp.UncNumber');
                narray = NET.GenericClass('Metas.UncLib.Core.Ndims.RealNArray', 'Metas.UncLib.LinProp.UncNumber');
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.LinAlg', {lu_res, narray, 'Metas.UncLib.LinProp.UncNumber'});
            end            
        end
        function l = LinAlg2(complex)
            if complex
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexLinAlg', {'Metas.UncLib.LinProp.UncNumber'});
            else
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.RealLinAlg', {'Metas.UncLib.LinProp.UncNumber'});
            end            
        end
        function l = NumLib(complex)
            if complex
                narray = NET.GenericClass('Metas.UncLib.Core.Ndims.ComplexNArray', 'Metas.UncLib.LinProp.UncNumber');
                number = NET.GenericClass('Metas.UncLib.Core.Complex', 'Metas.UncLib.LinProp.UncNumber');
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.NumLib', {narray, number});
            else
                narray = NET.GenericClass('Metas.UncLib.Core.Ndims.RealNArray', 'Metas.UncLib.LinProp.UncNumber');
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.NumLib', {narray, 'Metas.UncLib.LinProp.UncNumber'});
            end            
        end
        function l = NumLib2(complex)
            if complex
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNumLib', {'Metas.UncLib.LinProp.UncNumber'});
            else
                l = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNumLib', {'Metas.UncLib.LinProp.UncNumber'});
            end
        end
        function c = Double2ComplexNumber(d)
            c = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.Core.Number'});
            c.InitDblReIm(real(d), imag(d));
        end        
        function a = Double2Array(d)
            s = size(d);
            s = int32(s(:));
            if numel(d) == 0
                a = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
                a.InitNd(s);
            else
                d2 = d(:);
                if ~isreal(d)
                    a = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.Core.Number'});
                    a.InitDblReIm(real(d2), imag(d2));
                else
                    a = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.Core.Number'});
                    a.InitDbl(real(d2));
                end
                a.Reshape(s)
            end
        end
        function d = Convert2Double(x)
            if LinProp.IsArrayNet(x)
                s = int32(x.size);
                if LinProp.IsComplexNet(x)
                    d = complex(double(x.DblRealValue()), double(x.DblImagValue()));
                else
                    d = double(x.DblValue());
                end
                d = reshape(d, s);
            else
                if LinProp.IsComplexNet(x)
                    d = complex(x.DblRealValue(), x.DblImagValue());
                else
                    d = x.Value;
                end
            end
        end
        function m = Convert2UncArray(x)
            if x.IsArray
                m = x.NetObject;
            else
                if x.IsComplex
                    m = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                else
                    m = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});                    
                end
                m.Init2d(1, 1);
                m.SetItem2d(0, 0, x.NetObject);
            end 
        end
        function u = Convert2LinProp(x)
            if LinProp.IsArrayNet(x)
                if x.numel == 1
                    u = LinProp(x.GetItem2d(0, 0));
                else
                    u = LinProp(x);
                end
            else
                u = LinProp(x);
            end
        end
        function m = IndexMatrix(subs)
            n2 = numel(subs);
            s = zeros(1, n2);
            for i2 = 1:n2
                s(i2) = numel(subs{i2});
            end
            n1 = prod(s);
            m = zeros(n1, n2);
            temp = 1;
            for i2 = 1:n2
                temp_index = mod(floor((0:n1-1)./temp), s(i2)) + 1;
                m(:,i2) = subs{i2}(temp_index);
                temp = temp*s(i2); 
            end
        end
        function TF = IsComplexNet(x)
            TF = isa(x, 'Metas.UncLib.Core.Ndims.ComplexNArray<Metas*UncLib*LinProp*UncNumber>') || ...
                 isa(x, 'Metas.UncLib.Core.Complex<Metas*UncLib*LinProp*UncNumber>') || ...
                 isa(x, 'Metas.UncLib.Core.Ndims.ComplexNArray<Metas*UncLib*Core*Number>') || ...
                 isa(x, 'Metas.UncLib.Core.Complex<Metas*UncLib*Core*Number>');
        end
        function TF = IsArrayNet(x)
            TF = isa(x, 'Metas.UncLib.Core.Ndims.ComplexNArray<Metas*UncLib*LinProp*UncNumber>') || ...
                 isa(x, 'Metas.UncLib.Core.Ndims.RealNArray<Metas*UncLib*LinProp*UncNumber>') || ...
                 isa(x, 'Metas.UncLib.Core.Ndims.ComplexNArray<Metas*UncLib*Core*Number>') || ...
                 isa(x, 'Metas.UncLib.Core.Ndims.RealNArray<Metas*UncLib*Core*Number>');
        end
        function obj = XmlString2LinProp(s)
            UncPropLoadNETAssemblies('LinProp');
            try
                x = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                a = x.XmlDeserializeFromString(s);
            catch
                try
                    x = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    a = x.XmlDeserializeFromString(s);
                catch
                    try
                        x = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                        a = x.XmlDeserializeFromString(s);
                    catch
                        try
                            x = Metas.UncLib.LinProp.UncNumber();
                            a = x.XmlDeserializeFromString(s);
                        catch
                            error('Wrong structure of xml string')
                        end
                    end
                end
            end
            obj = LinProp(a);
        end
        function obj = XmlFile2LinProp(fp)
            UncPropLoadNETAssemblies('LinProp');
            try
                x = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                a = x.XmlDeserialize(fp);
            catch
                try
                    x = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    a = x.XmlDeserialize(fp);
                catch
                    try
                        x = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                        a = x.XmlDeserialize(fp);
                    catch
                        try
                            x = Metas.UncLib.LinProp.UncNumber();
                            a = x.XmlDeserialize(fp);
                        catch
                            error('Wrong structure of xml file')
                        end
                    end
                end
            end
            obj = LinProp(a);
        end
        function obj = BinaryFile2LinProp(fp)
            UncPropLoadNETAssemblies('LinProp');
            try
                x = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
                a = x.BinaryDeserialize(fp);
            catch
                try
                    x = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
                    a = x.BinaryDeserialize(fp);
                catch
                    try
                        x = NET.createGeneric('Metas.UncLib.Core.Complex', {'Metas.UncLib.LinProp.UncNumber'});
                        a = x.BinaryDeserialize(fp);
                    catch
                        try
                            x = Metas.UncLib.LinProp.UncNumber();
                            a = x.BinaryDeserialize(fp);
                        catch
                            error('Wrong structure of binary file')
                        end
                    end
                end
            end
            obj = LinProp(a);
        end
        function obj = System2LinProp(value, sys_inputs, sys_sensitivities)
            % Workaround to pass a SAFEARRAY with only one element
            if numel(sys_sensitivities) == 1
                sys_sensitivities = [sys_sensitivities 0];
            end
            sys_inputs = LinProp(sys_inputs);
            sys_inputs = LinProp.Convert2UncArray(sys_inputs);
            unc_number = Metas.UncLib.LinProp.UncNumber();
            unc_number.Init(value, sys_inputs.data, sys_sensitivities(:));
            obj = LinProp(unc_number);
        end
    end
end

function dispAsPages(name, value, isLoose)
    size_all = size(value);
    size_residual = size_all(3:end);
    size_residual = [size_residual ones(1, max(0, 2-length(size_residual)))];
    page_subscripts = cell(1, numel(size_residual));
    page_name = name;
    nPages = prod(size_residual);
    isComplex = ~isreal(value);
    for ii = 1:nPages
        [page_subscripts{:}] = ind2sub(size_residual,ii);
        page_values = value(:, :, page_subscripts{:});
        % Subscript assignment here removes the imag part if it is zero, so
        % we need to fix that for the call to disp.
        if isComplex
            page_values = complex(page_values);
        end

        if ~isempty(size_all(3:end))
            page_name = sprintf('%s(:,:,%s)', name, strjoin(strsplit(num2str(cell2mat(page_subscripts))), ','));
        end
        

        if (isLoose && ii==1); disp(' '); end
        disp([page_name ' = ']);
        if (isLoose); disp(' '); end
        disp(page_values);
        if (isLoose && ii ~= nPages); disp(' '); end
    end
end
