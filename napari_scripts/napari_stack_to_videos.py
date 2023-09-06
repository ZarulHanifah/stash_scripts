
## The contrast is actually inconsistent between stacks

import napari
import imageio
import skimage.io
import io

def export_as_mp4(viewer, layer_index, output_path, fps=30):
    layer = viewer.layers[layer_index]
    num_frames = layer.data.shape[0]
    original_height, original_width = layer.data.shape[1], layer.data.shape[2]
    
    # Adjust dimensions to be divisible by macro_block_size (e.g., 16)
    new_height = (original_height // 16) * 16
    new_width = (original_width // 16) * 16

    frames = []

    for frame_idx in range(num_frames):
        layer.visible = False  # Hide other layers
        layer.selected = True   # Select the current layer
        
        # Capture the current frame as a grayscale image
        current_frame = layer.data[frame_idx]
        
        # Normalize the frame to [0, 255] range and convert to uint8
        normalized_frame = ((current_frame - current_frame.min()) / (current_frame.max() - current_frame.min()) * 255).astype('uint8')
        
        # Resize the frame to match the desired dimensions
        resized_frame = normalized_frame[:new_height, :new_width]
        
        ### Save the frame using scikit-image's imsave
        buffer = io.BytesIO()
        skimage.io.imsave(buffer, resized_frame, format='png')
        buffer.seek(0)
        frame = imageio.imread(buffer)
        frames.append(frame)

        print(f"Exported frame {frame_idx + 1}/{num_frames}")

    # Export the frames as an MP4 video
    imageio.mimsave(output_path, frames, fps=fps)
    print(f"MP4 video saved to {output_path}")

layer_index = 0
output_path = 'output_video.mp4'
export_as_mp4(viewer, layer_index, output_path, fps=30)
