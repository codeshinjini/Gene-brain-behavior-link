function viewdicom(I, rot)
% Function to view dicom file with a scrollbar

if nargin < 1
    handles.rot = 0;
else
    handles.rot = rot;
end

handles.figImg = figure;
handles.imgI = imagesc(imrotate(I(:,:,1), rot), [min(I(:)) max(I(:))])
colorbar
%axis square
colormap jet

sizefig = get(handles.figImg, 'Position');

    handles.sld = uicontrol('Style', 'slider',...
            'Min',1,'Max',size(I,3),'Value',1,...
            'Position', [0 0 120 20],'SliderStep', [1/(size(I,3)-1) 10/(size(I,3)-1)],...
            'Callback', {@surfzlim,handles, I});
end

function surfzlim(hObj,event,handles, I)
    % Called to set zlim of surface in figure axes
    % when user moves the slider control
    slice = round(get(hObj,'Value'))
    set(handles.imgI, 'CData', imrotate(I(:,:,slice),handles.rot));

end