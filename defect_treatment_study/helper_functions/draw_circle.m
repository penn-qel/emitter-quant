function draw_circle(h, radius)
    global c;
    delete(c);
    c = viscircles([h(1), h(2)], radius);
end

