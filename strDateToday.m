
function stDate = strDateToday()

cc = clock;
stDate = [num2str(cc(1),'%02d') num2str(cc(2),'%02d') num2str(cc(3),'%02d')];
stDate = stDate(3:end);

end

