function [elf] = bcond1(elf)
    global im imm1 jm jmm1 fsm;

    elf(1, :) = elf(2, :); elf(im, :) = elf(imm1, :);
	elf(:, 1) = elf(:, 2); elf(:, jm) = elf(:, jmm1);
    elf = elf .* fsm;    
    return	
end