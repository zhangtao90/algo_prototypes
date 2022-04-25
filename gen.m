charv =  [ 'p' , 'v' , 'r' , 'a' , 'w' ];

gen_cmd1 = "syms ";

for i = 1:5
    for j = 1:5
        gen_cmd1 = gen_cmd1 + charv(i) + charv(j) + " ";
    end
end

gen_cmd2 = "Pk = [";

for i = 1:5
    for j = 1:5
        gen_cmd2 = gen_cmd2 + charv(i) + charv(j) + " ";
    end
    if(i==5)
        gen_cmd2 = gen_cmd2 + "];";
    else
        gen_cmd2 = gen_cmd2 + ";";
    end
end