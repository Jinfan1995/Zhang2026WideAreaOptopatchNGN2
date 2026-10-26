function [denoised_movie,U,S,V]=mini_ali_SVD_denoise(data,ncomp)
% mini_ali_SVD_denoise, compute sigular value decomposition (SVD) with 'ncomp'
% singular values
% input--
%   data: nrow x ncol x nframe image stacks
%   ncomp: number of kept componets
%   
% output--
%   denoised movie: reconsituted movie by U*S*V'
%   
% Created 04/21/2026 Jinfan Zhang
sz = size(data,[1:3]);

data = double(reshape(data,[],sz(3)));
[U,S,V] = svds(data,ncomp);
denoised_movie = U*S*V';
denoised_movie = single(reshape(denoised_movie,sz(1),sz(2),[]));
    
    