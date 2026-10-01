# PI Calculator

A super fast Python script I made that calculates digits of Pi (π) using binary splitting and has a cool telemetry code map like in satellite decoders.

Made by **[Dreamlexxi](https://github.com/Dreamlexxi)**!

## Why I made this

I wanted to see how fast I could calculate Pi in pure Python without having to install 50 different bloated libraries. At first, my simple loop took 10 seconds just to do 20,000 digits. So I rewrote it using Chudnovsky binary splitting with Python's big integers and `math.isqrt`, and now it can do **100,000 digits in about 0.3 seconds**.

I also added a neat "code map" grid in the terminal that lights up green as chunks finish decoding, kinda like satellite telemetry software.

## Features

- **Crazy fast:** 10,000 digits is practically instant (< 0.01s), 100k digits in ~0.3s, and 250k in ~1.4s.
- **Orbital Code Map:** A grid in your terminal that turns green (`■`) when chunks are done, yellow while working, and gray while waiting.
- **Clean Menu:** Everything stays on one screen instead of spamming your whole terminal scrollback.
- **Re-run button:** Just hit Enter or `R` to run it again with no hassle.
- **Find stuff in Pi:** Search for your birthday or any number inside Pi to see what decimal place it's at.
- **Zero pip installs:** Pure standard library. If you have Python 3.8+, you're good to go.

## How to run it

### Interactive Menu
Just run the script with nothing else:
```bash
python calculate_pi.py
```
This pops up the menu where you can pick presets (1k, 10k, 50k, 100k, 250k), switch algorithms, test your CPU, or search for numbers.

### CLI Mode
If you just want numbers fast or want to pipe it to a file:
```bash
# Get 10,000 digits
python calculate_pi.py 10000

# Get 100,000 digits and save it to a txt file
python calculate_pi.py 100000 -o pi.txt

# Watch the code map animate slowly (cool for demos)
python calculate_pi.py 500 --animate
```

## The Algorithms

I put 3 different methods in here so you can mess around with them:

1. **Chudnovsky (Binary Splitting)** - The default one. It's the same math supercomputers use for world records. It gets ~14 digits per step and does all the heavy lifting with pure integers.
2. **Gauss-Legendre** - Doubles the number of correct digits every loop iteration. It uses Python's `decimal` module. Really cool math, but a bit slower than pure integers once you pass 20k digits.
3. **Gregory-Leibniz** - The classic `4 * (1 - 1/3 + 1/5 - 1/7...)` series from math class. It's ridiculously slow (takes 100,000 steps just for 5 digits), but fun to benchmark against the big boys.

## Limits & Heads up (Don't crash your PC!)

Here are the actual limits and what happens if you push it:

- **Python's 4,300 digit limit:** Python 3.11 added a limit where it refuses to turn integers with more than 4,300 digits into strings. The script disables this with `sys.set_int_max_str_digits(0)` automatically, so you won't get that crash.
- **RAM usage:**
  - Up to **100,000 digits**: Uses almost nothing (~35 MB RAM).
  - **500,000 digits**: Around ~250 MB.
  - **1,000,000 digits**: Takes ~600 MB RAM and about 9-10 seconds.
  - **5,000,000+ digits**: You'll need at least 4 GB+ of free RAM. If you don't have enough RAM, Windows will start using swap memory and it'll slow down a lot.
- **Stack overflow?** Nope. The binary splitting tree splits in half every time, so even for 1,000,000 digits the recursion is only like 17 levels deep (Python's limit is 1,000).
- **Console lag:** If you calculate 100k digits, the script only prints a preview on screen so your Command Prompt doesn't freeze up trying to render 100,000 characters at once. Use `-o filename.txt` to get the full raw text!

## Benchmarks

On my machine, here is how fast it runs:

| Digits | Time |
|---|---|
| **1,000** | 0.0005s (instant) |
| **10,000** | 0.009s |
| **50,000** | 0.12s |
| **100,000** | 0.39s |
| **250,000** | 1.38s |

## License

MIT License. Feel free to copy, modify, or do whatever you want with it!

Made by [Dreamlexxi](https://github.com/Dreamlexxi) on GitHub.
