function x = morse_decode(sir)
  m = morse();
  for k = 1:length(sir)
    if sir(k) == "."
      m = m{2};
    elseif sir(k) == "-"
      m = m{3};
    end
    if isempty(m)
      ch = "*";
      return
    end
  end
  x = m{1};
end