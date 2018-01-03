function x = multiple_decode(sir)
    % decode('string of dots, dashes and spaces')
    x = [];
    sir = [sir ' '];
    while ~isempty(sir);
       k = find(sir == ' ',1);
       ch = morse_decode(sir(1:k));
       x = [x ch];
       code(1:k) = [];
       % Many blanks in the code is worth one in the text.
       if ~isempty(code) && code(1) == ' '
           text = [x ' '];
          while ~isempty(sir) && sir(1) == ' '
             sir(1) = [];
          end
       end
    end 
  
