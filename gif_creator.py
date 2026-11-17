import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import tqdm

L = 128
B = [0.600, 0.700, 0.441, 0.442, 0.399, 0.440, 0.443, 0.000, 0.200, 0.300, 1.000]

for valor in B:
    filename = f"results/spins_gif_L{L:03}_B{valor:.3f}.bin"

    with open(filename, "rb") as f:
        data = np.fromfile(f, dtype=np.int8)

    total_frames = len(data) // (L * L)
    complete_frames = data[:total_frames * L**2].reshape(-1, L, L)

    # Keep early thermalization frames dense, then subsample later frames
    culling_threshold = int(L / 4)
    frames = np.concatenate([
        complete_frames[:culling_threshold],
        complete_frames[culling_threshold::10]
    ])

    fig = plt.figure(figsize=(5.12, 5.12), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    im = ax.imshow(frames[0], cmap='binary', interpolation='nearest', animated=True)

    def update(frame_idx):
        im.set_array(frames[frame_idx])
        return [im]

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=tqdm.tqdm(range(len(frames)), desc=f"Rendering B={valor:.3f}"),
        interval=100,
        blit=True,
    )

    ani.save(f"plots/ising_thermalization_L{L:03}_B{valor:.3f}_final_slowed_at_first_binary.gif", writer='pillow')
    plt.close()
