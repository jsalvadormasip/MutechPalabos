import os
from PIL import Image

def clear_folder(folder_path):
    """
    Delete all files in the specified folder.
    """
    if os.path.exists(folder_path):
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if os.path.isfile(file_path):
                os.unlink(file_path)

def convert_ppm_to_gif(input_folder, output_folder):
    """
    Convert all PPM files in the input folder to GIF files and save them to the output folder.
    """
    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # List all files in the input folder
    for filename in os.listdir(input_folder):
        # Check if the file is a PPM file
        if filename.lower().endswith('.ppm') and filename.lower().startswith('vnorm_x'):
            # Construct full file path
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename[:-4] + '.gif')
            
            # Open the PPM image
            with Image.open(input_path) as img:
                # Save the image as a GIF
                img.save(output_path, 'GIF')
            
            print(f"Converted {input_path} to {output_path}")

def create_animated_gif_from_folder(source_folder, output_folder, output_filename, duration=500, loop=0):
    """
    Create an animated GIF from all GIF files in a source folder and save it to an output folder.
    """
    # Get a list of all GIF files in the source folder
    gif_files = [os.path.join(source_folder, f) for f in os.listdir(source_folder) if f.lower().endswith('.gif')]
    
    # Sort the files to ensure proper order
    gif_files.sort()
    
    # Ensure there are GIF files to process
    if not gif_files:
        print("No GIF files found in the source folder.")
        return

    # Open the first image and convert it to RGBA mode
    frames = [Image.open(gif_files[0]).convert('RGBA')]

    # Append the rest of the images
    for gif_file in gif_files[1:]:
        img = Image.open(gif_file).convert('RGBA')
        frames.append(img)

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Save the frames as an animated GIF
    output_path = os.path.join(output_folder, output_filename)
    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=loop,
        disposal=2  # To properly handle transparency
    )

    print(f"Animated GIF saved to {output_path}")

if __name__ == "__main__":
    input_folder = 'C:/Users/jordi/Desktop/MutechInternship/Palabos/palabos-master/examples/showCases/jordiPowerFlowCopy/tmp'  # Replace with the path to your folder containing PPM files
    intermediate_gif_folder = 'C:/Users/jordi/Desktop/MutechInternship/gifs'  # Temporary folder for converted GIFs
    output_folder = 'C:/Users/jordi/Desktop/MutechInternship/AnimatedGif'  # Folder for the final animated GIF
    output_filename = 'animated.gif'  # Name of the final animated GIF

    # Clear the GIF and animated GIF folders
    clear_folder(intermediate_gif_folder)
    clear_folder(output_folder)

    # Step 1: Convert all PPM files to GIF files
    convert_ppm_to_gif(input_folder, intermediate_gif_folder)

    # Step 2: Create an animated GIF from the converted GIF files
    create_animated_gif_from_folder(intermediate_gif_folder, output_folder, output_filename)
