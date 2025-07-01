using Plots 

gr()  # Oder plotly(), pyplot(), etc.

# Erzeuge 4 Subplots
p1 = plot(sin, 0, 2π, title="Sinus", label="sin(x)")
p2 = plot(cos, 0, 2π, title="Cosinus", label="cos(x)")
p3 = plot(tan, 0, π/2 - 0.1, title="Tangens", label="tan(x)")
p4 = plot(x -> x^2, 0, 3, title="x^2", label="x²")

# Kombiniere sie zu einem 2×2-Raster
display(plot(p1, p2, p3, p4, layout=(2,2), size=(800,600), legend=:bottomright))

# Erzeuge 4 Subplots
p5 = plot(sin, 0, 2π, title="Sinus", label="sin(x)")
p6 = plot(cos, 0, 2π, title="Cosinus", label="cos(x)")
p7 = plot(tan, 0, π/2 - 0.1, title="Tangens", label="tan(x)")
p8 = plot(x -> x^2, 0, 3, title="x^2", label="x²")

# Kombiniere sie zu einem 2×2-Raster
display(plot(p5, p6, p7, p8, layout=(2,2), size=(800,600), legend=:bottomright))