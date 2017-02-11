function cor=init_coriolis(grid)
    global im jm kb
    cor=create_field(zeros(im,jm,kb), grid, 3);
end