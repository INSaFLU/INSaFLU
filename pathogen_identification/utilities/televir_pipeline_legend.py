import matplotlib.pyplot as plt

pipeline_steps_to_r_colours = {
    "root": "lightblue",
    "Extra QC": "cadetblue",
    "Viral Enrichment": "darkgreen",
    "Host Depletion": "blueviolet",
    "Assembly": "brown",
    "Contig Classification": "darkorange",
    "Read Classification": "deeppink",
    "Metagenomics Combine": "dodgerblue",
    "Remapping": "khaki",
    "Remap Filtering": "darkslategray",
    "Map Filtering": "darkolivegreen",
    "Request Mapping": "darkolivegreen",
    "Leaves": "lightblue",
}

# Reverse the order of the dictionary
reversed_dict = dict(reversed(pipeline_steps_to_r_colours.items()))

# Create a figure and axis
fig, ax = plt.subplots()

# Set the aspect ratio to equal
ax.set_aspect("equal")

# Set the spacing between circles
spacing = 1.5

# Iterate over the reversed dictionary items
for i, (desc, color) in enumerate(reversed_dict.items()):
    # Calculate the y-coordinate with spacing
    y = i * spacing

    # Draw a circle with the specified color
    circle = plt.Circle((0, y), 0.5, color=color, fill=True)

    # Add the circle to the plot
    ax.add_artist(circle)

    # Add the description next to the circle
    ax.text(1, y, desc, va="center")

# Set the x and y limits
ax.set_xlim(-1, 11)
ax.set_ylim(-1, len(reversed_dict) * spacing)

# Remove the axis ticks
ax.set_xticks([])
ax.set_yticks([])

# Save the plot as a JPG image with increased resolution
fig.savefig("pipeline_steps.jpg", bbox_inches="tight", dpi=500)

# Show the plot
plt.show()
