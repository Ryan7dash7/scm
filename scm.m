function chi = scm(t,p,w,tau_min,tau_max)
%SCM Spatial correlation metric.
%   CHI = SCM(T,P,W,TAU_MIN,TAU_MAX) computes the spatial correlation
%   metric of pressure P at time T as described in "Analysis of Axial
%   Compressor Stall Inception Using Unsteady Casing Pressure Measurements"
%   by Cameron and Morris. The spatial correlation metric is based on a
%   windowed, two-point correlation function between adjacent sensors and
%   provides a scalar function that is nonzero only when disturbances that
%   rotate around the compressor annulus in the direction of the rotor's
%   rotation are present.
%
%   W, the cross-correlation window length, specifies the time scale
%   (length scale for rotating disturbances) at which the "comparison" of
%   data traces is made. For example, W = 1/N will cause CHI to be
%   sensitive to rotating disturbances with length scale on the order of
%   the compressor circumference and smaller. TAU_MIN and TAU_MAX, the time
%   lags, specify the desired rotational speed range. Disturbances that are
%   rotating outside of the specified range of angular velocity will not
%   contribute to the net CHI value.
%
%   For column-oriented data analysis, use SCM(T',P',W,TAU_MIN,TAU_MAX).

%   Version 2017.03.14
%   Copyright 2017 Ryan McGowan

[n_theta,n_t] = size(p);
chi = zeros(n_theta,n_t);
h = waitbar(0,'Please wait...');
for i = 1:n_theta
    for j = 1:n_t
        % Window home trace around t_0
        t_0 = t(j);
        W = heaviside(t-(t_0-w/2))-heaviside(t-(t_0+w/2));
        pHat = p(i,:).*W;
        
        % Compute cross-correlation between windowed signal and neighbor
        % time trace
        if i == n_theta
            [R,lags] = xcorr(pHat,p(1,:));
        else
            [R,lags] = xcorr(pHat,p(i+1,:));
        end
        
        % Compute SCM
        tau = lags*mean(diff(t));
        [~,i_tau_min] = min(abs(tau-tau_min));
        [~,i_tau_max] = min(abs(tau-tau_max));
        i_tau = i_tau_min:i_tau_max;
        chi(i,j) = trapz(tau(i_tau),R(i_tau));
        
        waitbar(((i-1)*n_t+j)/(n_theta*n_t),h)
    end
end
close(h)

end

