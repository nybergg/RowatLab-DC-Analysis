function Image = RasterizeText(String)

% Image = RasterizeText(String,Font,FontSize)
%
%   Creates a monochrome image with a rasterized version of the text
%   specified in String with font specified by Font. Font is a structure 
%   previously produced by BitmapFont.m which is saved in the current 
%   directory.
%
%   No line wrapping occurs, but the function can process the newline
%   character. Output size is unpredictable without first analysing the
%   font. Best results are likely with fixed width fonts.
%
% Daniel Warren
% Particle Therapy Cancer Research Institute
% University of Oxford

if ~exist('String','var')
    String = 'No string specified.';
end

% Preprocess text. Only allowing two types of whitespace: \n and space
% Replace tab with four spaces. Remove all other ASCII control characters.
String = strrep(String,sprintf('\t'),sprintf('    '));
String = strrep(String,sprintf('\r\n'),sprintf('\n'));
ControlChars = sprintf('%c',[0:9 11:31 127]);
for i = 1:length(ControlChars)
    String(String==ControlChars(i)) = [];
end

% Create a rasterized font
Characters = unique(String(String ~= ' ' & String ~= sprintf('\n')));
Font = importdata('imageFonts.mat');

Image = logical([]); % This array will grow as the output is built up

l = 0; % Line number - starts at 0
x = 0; % X location - starts a 0, but 0 will never be written to

SpaceSize = ceil(0.33*Font.Size);

for i = 1:length(String)
    switch String(i)
        case ' '
            % Avoid overwriting parts of characters below the baseline on
            % the line above by only assigning one element.
            Image(l*Font.Size + Font.Size, x + SpaceSize) = false;
            x = x+SpaceSize;
        case sprintf('\n')
            l = l+1;
            % Unnecessary to grow array, but could help speed. Again, only
            % assign one element.
            Image(l*Font.Size + Font.Size, size(Image,2)) = false;
            x = 0;
        otherwise
            index = Font.Characters==String(i);
            CharSize = size(Font.Bitmaps{index});
            % Grow array so can perform boolean OR, which will avoid
            % background of character overwriting characters extending
            % below the baseline on the line above.
            Image(l*Font.Size + CharSize(1), x + CharSize(2)) = false;
            Image(l*Font.Size + (1:CharSize(1)), x + (1:CharSize(2))) = ...
                Image(l*Font.Size + (1:CharSize(1)), x + (1:CharSize(2))) | Font.Bitmaps{index};
            x = x+CharSize(2);
    end
end

end