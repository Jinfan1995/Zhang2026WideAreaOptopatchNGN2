function [shift,dfdx,dfdy,reference_image]=mini_ali_motion_corr_info(movie, reference_frame)

reference_image = mean(single(movie(:,:,reference_frame)), 3);

dfdx = (reference_image(:,3:end) - reference_image(:,1:end-2)) / 2;
dfdy = (reference_image(3:end,:) - reference_image(1:end-2,:)) / 2;

dfdx = [dfdx(:,1), dfdx, dfdx(:,end)];
dfdy = [dfdy(1,:); dfdy; dfdy(end,:)];


ss = size(movie);
shift = zeros(2, ss(3), 'single');

chunkSize = 200;

A = single([dfdx(:), dfdy(:)]);

AtA = A' * A;

At = A';

for t0 = 1:chunkSize:ss(3)
    t1 = min(t0 + chunkSize - 1, ss(3));
    frames = reshape(single(movie(:,:,t0:t1)), [], t1-t0+1);
    shift_chunk = AtA \ (At * frames);
    shift(:,t0:t1) = shift_chunk;
end

