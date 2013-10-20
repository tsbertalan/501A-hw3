function [] = bump_keypress_callback()

c = get(gcf, 'CurrentCharacter');

set(gcf, 'UserData', c);

