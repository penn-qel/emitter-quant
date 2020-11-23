function [hsv_img, hue, sat, val] = calc_hsv_from_fft(im_fft)
%CALC_HSV Calculates the hsv values from the fft
fft1 = im_fft{1};
fft2 = im_fft{2};

hue = mod(angle(fft2), pi) / pi; %Angle

sat = abs(fft2) ./ (abs(fft1) + abs(fft2)); %Visibility
sat = 2 * sat; %NB: From half angle formula, max for dipole is s = 0.5
if max(sat(:) > 1) %In case numerical issues make s larger than 1 for some cases
    warning('Warning - detected sat > 1, cutting this off');
    warning(['Max sat: ' num2str(max(sat(:)))]);
    sat = min(sat,1);
end

val = abs(fft1) / max(abs(fft1(:))); %Average value, looks better
%{
v_ref = max(abs(fft1(:)));
val = 2*(1 - v_ref./(v_ref + abs(fft1))); %Stereographic, may have more information
%}

hsv_img = cat(3, hue, sat, val);

end

