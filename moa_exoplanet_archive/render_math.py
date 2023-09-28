import re
from pathlib import Path

import pylab


def render_math(math_string: str, entry_number: int):
    math_renders_root_directory = Path('math_renders')
    math_renders_png_directory = math_renders_root_directory.joinpath('png')
    math_renders_svg_directory = math_renders_root_directory.joinpath('svg')
    formula = rf'${math_string}$'

    fig = pylab.figure()
    text = fig.text(0, 0, formula)

    # Saving the figure will render the text.
    dpi = 300
    fig.savefig('temporary.png', dpi=dpi)

    # Now we can work with text's bounding box.
    bbox = text.get_window_extent()
    width, height = bbox.size / float(dpi) + 0.005
    # Adjust the figure size so it can hold the entire text.
    fig.set_size_inches((width, height))

    # Adjust text's vertical position.
    dy = (bbox.ymin / float(dpi)) / height
    text.set_position((0, -dy))

    # Save the adjusted text.
    math_renders_png_directory.mkdir(exist_ok=True, parents=True)
    math_renders_svg_directory.mkdir(exist_ok=True, parents=True)
    file_name_string = (math_string.replace('\\', 'backslash')
                        .replace('/', 'slash').replace('%', 'percent')
                        .replace('{', 'leftbrace').replace('}', 'rightbrace')
                        .replace('.', 'dot').replace('>', 'gt').replace('<', 'lt')
                        .replace('+', 'plus').replace('-', 'minus')
                        .replace('^', 'caret').replace(',', 'comma')
                        .replace(' ', '_').replace('(', 'leftparen')
                        .replace(')', 'rightparen').replace('=', 'equals')
                        )
    math_render_png_path = math_renders_png_directory.joinpath(f'{entry_number:02}_{file_name_string}.png')
    fig.savefig(math_render_png_path, dpi=dpi)
    math_render_svg_path = math_renders_svg_directory.joinpath(f'{entry_number:02}_{file_name_string}.svg')
    fig.savefig(math_render_svg_path, dpi=dpi)
    pylab.close(fig)


def render_all_math():
    documentation_path = Path('ipac_documentation/documentation.md')
    with documentation_path.open() as documentation_file:
        documentation_content = documentation_file.read()
        search_pattern = rf'\$([^\n\$]+)\$'
        match_list = re.findall(search_pattern, documentation_content)
        for match_index, match in enumerate(match_list):
            render_math(match, match_index)


if __name__ == '__main__':
    render_all_math()
