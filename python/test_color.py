for i in range(256):
    print(f"\x1B[38;5;{i}m {i:3d} \033[0m", end=' ')
    if i % 16 == 15:
        print()

for _ in range(10) : print()

r = 255
g = 0
b = 0
print(f"\x1B[38;2;{r};{g};{b}m ({r:3d},{g:3d},{b:3d})\033[0m", end=' ')
