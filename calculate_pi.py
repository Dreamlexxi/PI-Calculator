#!/usr/bin/env python3
"""
calculate_pi.py

A high-performance arbitrary-precision Pi (π) calculator using pure-integer
Binary Splitting Chudnovsky, dynamic-precision Gauss-Legendre, and paired
Gregory-Leibniz algorithms.

Author: Dreamlexxi (https://github.com/Dreamlexxi)
License: MIT
"""

import argparse
import ctypes
import math
import os
import sys
import time
from decimal import Decimal, getcontext

# Lift Python 3.11+ integer-to-string digit conversion ceiling
if hasattr(sys, "set_int_max_str_digits"):
    try:
        sys.set_int_max_str_digits(0)
    except Exception:
        pass

# Ensure console supports UTF-8 characters without encoding crashes
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# Enable virtual terminal ANSI sequence processing on Windows consoles
if os.name == "nt":
    try:
        kernel32 = ctypes.windll.kernel32
        kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
    except Exception:
        pass


def clear_screen():
    """Clear the console screen for single-page view transitions."""
    if sys.stdout.isatty():
        os.system("cls" if os.name == "nt" else "clear")


class OrbitalCodeMap:
    """
    Real-time orbital telemetry code map visualization.
    Renders computing frames across a multi-channel matrix:
      ■ Green  : Frame Decoded & Verified [OK]
      ■ Yellow : Processing Active Subtree [*]
      · Gray   : Frame Pending Transmission [ ]
    """

    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    CYAN = "\033[96m"
    GRAY = "\033[90m"
    WHITE = "\033[97m"
    BOLD = "\033[1m"
    RESET = "\033[0m"

    def __init__(self, total_steps: int, total_blocks: int = 100, cols: int = 25):
        self.total_steps = max(1, total_steps)
        self.total_blocks = total_blocks
        self.cols = cols
        self.rows = (total_blocks + cols - 1) // cols
        self.lines_printed = 0
        self.is_tty = sys.stdout.isatty()
        self.start_time = time.perf_counter()
        self.last_update_time = 0.0

    def update(self, current_step: int, extra: str = "", force: bool = False):
        """Update telemetry frames and console metrics."""
        now = time.perf_counter()
        if not force and (now - self.last_update_time < 0.03) and current_step < self.total_steps:
            return

        self.last_update_time = now
        elapsed = now - self.start_time
        percent = min(100.0, (current_step / self.total_steps) * 100)
        active_block = int((percent / 100.0) * self.total_blocks)

        if percent > 0 and percent < 100:
            eta = (elapsed / (percent / 100.0)) - elapsed
            eta_str = f"{eta:.2f}s"
        elif percent >= 100:
            eta_str = "0.00s"
        else:
            eta_str = "--"

        lines = []
        lines.append(f"{self.CYAN}┌── [ ORBITAL CODE MAP - TELEMETRY DECODER ] ───────────────────┐{self.RESET}")

        for r in range(self.rows):
            row_cells = []
            for c in range(self.cols):
                b_idx = r * self.cols + c
                if b_idx < self.total_blocks:
                    if b_idx < active_block or percent >= 100.0:
                        row_cells.append(f"{self.GREEN}■{self.RESET}")
                    elif b_idx == active_block:
                        row_cells.append(f"{self.YELLOW}■{self.RESET}")
                    else:
                        row_cells.append(f"{self.GRAY}·{self.RESET}")
            lines.append(f"{self.CYAN}│{self.RESET}  " + " ".join(row_cells) + f"  {self.CYAN}│{self.RESET}")

        lines.append(f"{self.CYAN}└──────────────────────────────────────────────────────────────┘{self.RESET}")
        lines.append(
            f"  Frames: {self.GREEN}■{self.RESET} [OK] Synced  "
            f"{self.YELLOW}■{self.RESET} [*] In-Lock  "
            f"{self.GRAY}·{self.RESET} [ ] Pending"
        )
        lines.append(
            f"  {self.BOLD}Lock:{self.RESET} [{self.GREEN}{percent:5.1f}%{self.RESET}]  "
            f"Frame: {current_step}/{self.total_steps}  "
            f"Elapsed: {elapsed:.2f}s  ETA: {eta_str}  {extra}"
        )

        if self.is_tty:
            if self.lines_printed > 0:
                sys.stdout.write(f"\033[{self.lines_printed}F")
                for line in lines:
                    sys.stdout.write(f"\033[K{line}\n")
            else:
                for line in lines:
                    sys.stdout.write(f"{line}\n")
            sys.stdout.flush()
            self.lines_printed = len(lines)
        else:
            if force:
                print("\n".join(lines))

    def finish(self):
        """Lock telemetry frame grid upon completion."""
        self.update(self.total_steps, extra="[CARRIER LOCKED / 100%]", force=True)
        print()

def _bs_chudnovsky(a: int, b: int) -> tuple:
    """
    Evaluates P(a,b), Q(a,b), T(a,b) over integer range [a, b).
    Avoids floating-point math completely by computing large integers directly.
    """
    if b - a == 1:
        if a == 0:
            return 1, 1, 13591409
        P = (6 * a - 5) * (2 * a - 1) * (6 * a - 1)
        Q = a**3 * 10939058860032000
        T = P * (13591409 + 545140134 * a)
        return P, Q, -T if a % 2 == 1 else T

    m = (a + b) // 2
    P1, Q1, T1 = _bs_chudnovsky(a, m)
    P2, Q2, T2 = _bs_chudnovsky(m, b)
    return P1 * P2, Q1 * Q2, T1 * Q2 + P1 * T2


def _merge_binary_chunks(chunks: list) -> tuple:
    """Balanced tree reduction for chunked integer tuples."""
    while len(chunks) > 1:
        new_chunks = []
        for i in range(0, len(chunks), 2):
            if i + 1 < len(chunks):
                P1, Q1, T1 = chunks[i]
                P2, Q2, T2 = chunks[i + 1]
                new_chunks.append((P1 * P2, Q1 * Q2, T1 * Q2 + P1 * T2))
            else:
                new_chunks.append(chunks[i])
        chunks = new_chunks
    return chunks[0]


def calculate_pi_chudnovsky(
    digits: int, code_map: OrbitalCodeMap = None, step_delay: float = 0.0
) -> str:
    """
    Calculate Pi to `digits` places using integer binary splitting.
    """
    terms = int(digits / 14.181647462725477) + 1

    if code_map:
        code_map.update(0, extra="Aligning carrier sequence...")
        num_chunks = min(terms, 100)
        chunk_size = (terms + num_chunks - 1) // num_chunks
        chunk_results = []

        for i in range(num_chunks):
            a = i * chunk_size
            b = min(terms, (i + 1) * chunk_size)
            if a >= b:
                continue
            res = _bs_chudnovsky(a, b)
            chunk_results.append(res)

            elapsed = max(0.0001, time.perf_counter() - code_map.start_time)
            rate = int((b * 14.1816) / elapsed)
            code_map.update(i + 1, extra=f"Rate: {rate:,.0f} dig/s")
            if step_delay > 0:
                time.sleep(step_delay)

        P, Q, T = _merge_binary_chunks(chunk_results)
    else:
        # Maximum speed direct binary splitting (no chunking overhead)
        P, Q, T = _bs_chudnovsky(0, terms)

    # Scaled integer square root via C-accelerated math.isqrt
    extra = 10
    total_D = digits + extra
    sqrt_c = math.isqrt(10005 * 10**(2 * total_D))
    pi_scaled = (426880 * sqrt_c * Q) // T

    # Round to nearest target decimal place
    pi_rounded = (pi_scaled + 5 * 10**(extra - 1)) // (10**extra)
    pi_s = str(pi_rounded)
    res = pi_s[0] + "." + pi_s[1:digits + 1]

    if code_map:
        code_map.finish()

    return res

def calculate_pi_gauss_legendre(
    digits: int, code_map: OrbitalCodeMap = None, step_delay: float = 0.0
) -> str:
    """
    Calculate Pi using Gauss-Legendre (Salamin-Brent) algorithm.
    """
    iters = digits.bit_length() + 2

    if code_map:
        code_map.update(0, extra="Initiating AGM convergence...")

    getcontext().prec = digits + 10

    a = Decimal(1)
    b = Decimal(1) / Decimal(2).sqrt()
    t = Decimal(1) / Decimal(4)
    p = Decimal(1)

    for i in range(1, iters + 1):
        an = (a + b) / 2
        bn = (a * b).sqrt()
        t -= p * (a - an) ** 2
        a, b = an, bn
        p *= 2

        if code_map:
            code_map.update(i, extra=f"Iteration {i}/{iters}")
        if step_delay > 0:
            time.sleep(step_delay)

    pi_val = ((a + b) ** 2) / (4 * t)
    getcontext().prec = digits
    result = str(+pi_val)

    if code_map:
        code_map.finish()

    return result

def calculate_pi_leibniz(
    iterations: int = 100000, code_map: OrbitalCodeMap = None, step_delay: float = 0.0
) -> str:
    """
    Approximates Pi using paired-term Gregory-Leibniz summation.
    """
    half_iters = iterations // 2
    update_interval = max(1, half_iters // 100)
    acc = 0.0

    if code_map:
        code_map.update(0, extra="Decoding series pairs...")

    for k in range(half_iters):
        acc += 2.0 / ((4 * k + 1) * (4 * k + 3))

        if code_map and (k % update_interval == 0 or k == half_iters - 1):
            step = int((k / half_iters) * 100) + 1
            code_map.update(step, extra=f"Pairs {2*k+2:,}/{iterations:,}")
            if step_delay > 0:
                time.sleep(step_delay)

    pi_approx = 4.0 * acc

    if code_map:
        code_map.finish()

    return f"{pi_approx:.15f}"

def run_single_calculation(
    digits: int,
    algorithm: str = "chudnovsky",
    use_code_map: bool = True,
    delay: float = 0.0,
    output_file: str = None,
) -> tuple:
    """Run calculation, rendering orbital telemetry code map."""
    clear_screen()
    print("=" * 66)
    print("            PI (π) CALCULATOR - ORBITAL TELEMETRY DECODER         ")
    print("                         by Dreamlexxi                            ")
    print("=" * 66)
    print(f"Algorithm: {algorithm.capitalize()}")
    print(f"Target:    {digits:,} decimal digits")
    if delay > 0:
        print(f"Pacing:    {delay:.3f}s frame delay")
    print()

    code_map = None
    if use_code_map:
        if algorithm == "chudnovsky":
            terms = int(digits / 14.181647462725477) + 1
            total_steps = min(terms, 100)
        elif algorithm == "gauss-legendre":
            total_steps = digits.bit_length() + 2
        elif algorithm == "leibniz":
            total_steps = 100
        code_map = OrbitalCodeMap(total_steps=total_steps, total_blocks=100, cols=25)

    start_time = time.perf_counter()

    if algorithm == "chudnovsky":
        pi_str = calculate_pi_chudnovsky(digits, code_map=code_map, step_delay=delay)
    elif algorithm == "gauss-legendre":
        pi_str = calculate_pi_gauss_legendre(digits, code_map=code_map, step_delay=delay)
    elif algorithm == "leibniz":
        iterations = max(digits, 100000)
        pi_str = calculate_pi_leibniz(iterations, code_map=code_map, step_delay=delay)

    elapsed = time.perf_counter() - start_time

    print("Decoded Value:")
    print("-" * 66)
    if len(pi_str) > 1000:
        print(pi_str[:500] + "\n\n... [output truncated for terminal display] ...\n\n" + pi_str[-100:])
    else:
        print(pi_str)
    print("-" * 66)
    print(f"Computation completed in {elapsed:.4f} seconds.")

    if output_file:
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(pi_str)
        print(f"Saved full stream to: {output_file}")

    return pi_str, elapsed


def benchmark_algorithms(digits: int = 10000):
    """Side-by-side performance benchmark of all supported algorithms."""
    clear_screen()
    print("=" * 66)
    print(f"        ALGORITHM BENCHMARK & COMPARISON ({digits:,} digits)       ")
    print("=" * 66)
    print("Benchmarking unthrottled computational pipelines...\n")

    t0 = time.perf_counter()
    pi_chud = calculate_pi_chudnovsky(digits, code_map=None, step_delay=0.0)
    t_chud = time.perf_counter() - t0

    t0 = time.perf_counter()
    pi_gauss = calculate_pi_gauss_legendre(digits, code_map=None, step_delay=0.0)
    t_gauss = time.perf_counter() - t0

    t0 = time.perf_counter()
    pi_leib = calculate_pi_leibniz(iterations=1000000, code_map=None, step_delay=0.0)
    t_leib = time.perf_counter() - t0

    print("┌─────────────────┬──────────────┬──────────────┬──────────────────┐")
    print("│ Algorithm       │ Time (sec)   │ Precision    │ Convergence      │")
    print("├─────────────────┼──────────────┼──────────────┼──────────────────┤")
    print(f"│ Chudnovsky (BS) │ {t_chud:10.5f} s │ {digits:6d} dig  │ ~14.2 dig/term   │")
    print(f"│ Gauss-Legendre  │ {t_gauss:10.5f} s │ {digits:6d} dig  │ Quadratic (2x)   │")
    print(f"│ Gregory-Leibniz │ {t_leib:10.5f} s │ ~6 digits    │ Linear (O(1/n))  │")
    print("└─────────────────┴──────────────┴──────────────┴──────────────────┘")

    chud_str = pi_chud[:min(digits, len(pi_chud))]
    gauss_str = pi_gauss[:min(digits, len(pi_gauss))]
    if chud_str == gauss_str:
        print("\033[92m✔ Verification pass:\033[0m Chudnovsky and Gauss-Legendre match 100%!")
    else:
        print("\033[93m⚠ Divergence observed at trailing guard digits.\033[0m")
    print("-" * 66)


def search_in_pi(pi_str: str):
    """Search for a specific sequence (birthday, callsign, pattern) in Pi digits."""
    clear_screen()
    print("=" * 66)
    print("               PATTERN SEARCH IN PI (π) STREAM                   ")
    print("=" * 66)

    if not pi_str:
        print("No Pi digits loaded. Run a calculation first!")
        return

    print(f"Active buffer: {len(pi_str) - 2:,} decimal places.")
    query = input("\nEnter digit sequence to locate (e.g. 1415, 2026, MMDD): ").strip()
    if not query or not query.isdigit():
        print("Invalid sequence. Numbers only.")
        return

    start_offset = 2 if pi_str.startswith("3.") else 0
    pos = pi_str.find(query, start_offset)

    if pos != -1:
        decimal_idx = pos - start_offset + 1
        start_sub = max(start_offset, pos - 8)
        end_sub = min(len(pi_str), pos + len(query) + 8)
        context_str = (
            pi_str[start_sub:pos]
            + f"\033[92m[{query}]\033[0m"
            + pi_str[pos + len(query):end_sub]
        )
        print(f"\n\033[92m✔ Sequence located!\033[0m")
        print(f"  Decimal Position: Place #{decimal_idx:,}")
        print(f"  Context Frame:    ...{context_str}...")
    else:
        print(f"\n\033[93m✘ Sequence '{query}' not found\033[0m in the first {len(pi_str) - start_offset:,} digits.")
        print("  Increase target precision to search deeper into the sequence.")


def post_calc_actions(
    digits: int,
    algorithm: str,
    use_code_map: bool,
    delay: float,
    pi_str: str,
) -> str:
    """Action prompt after calculation completes."""
    while True:
        print()
        print("┌──────────────────────────────────────────────────────────────┐")
        print("│  [Enter / R] Re-run  │  [S] Save to File  │  [F] Search in π │")
        print("│  [M] Main Menu       │  [Q] Exit                             │")
        print("└──────────────────────────────────────────────────────────────┘")
        choice = input("Select [default: R]: ").strip().lower()

        if choice in ("", "r", "rerun"):
            return "rerun"
        elif choice in ("s", "save"):
            fname = input(f"Filename [default: pi_{digits}_digits.txt]: ").strip()
            if not fname:
                fname = f"pi_{digits}_digits.txt"
            with open(fname, "w", encoding="utf-8") as f:
                f.write(pi_str)
            print(f"\033[92m✔ Stream exported to {fname}\033[0m")
        elif choice in ("f", "find", "search"):
            search_in_pi(pi_str)
        elif choice in ("m", "menu"):
            return "menu"
        elif choice in ("q", "quit", "exit"):
            return "quit"
        else:
            print("Invalid input. Choose R, S, F, M, or Q.")


def interactive_menu():
    """Main terminal console dashboard."""
    digits = 10000
    algorithm = "chudnovsky"
    use_code_map = True
    delay = 0.0
    last_pi_str = ""
    status_msg = ""

    while True:
        clear_screen()
        anim_status = f"ON ({delay:.3f}s)" if delay > 0 else "OFF (Max Speed)"
        map_status = "ON" if use_code_map else "OFF"

        print("=" * 66)
        print("            PI (π) CALCULATOR - ORBITAL TELEMETRY DECODER         ")
        print("                    Created by Dreamlexxi on GitHub               ")
        print("=" * 66)
        if status_msg:
            print(f"  {status_msg}\n" + "-" * 66)
            status_msg = ""

        print(f"  [1] Decode Pi Stream         (Target: {digits:,} digits, {algorithm.capitalize()})")
        print(f"  [2] Target Precision         (Current: {digits:,})")
        print(f"  [3] Decoding Algorithm       (Current: {algorithm.capitalize()})")
        print(f"  [4] Code Map & Telemetry     (Map: {map_status}, Pacing: {anim_status})")
        print("  [5] Benchmark Algorithms")
        print("  [6] Search Pattern in Pi Stream")
        print("  [7] Save Stream to File")
        print("  [0] Exit")
        print("-" * 66)

        try:
            choice = input("Select an option [0-7, default: 1]: ").strip()
        except (KeyboardInterrupt, EOFError):
            clear_screen()
            print("Session terminated. Goodbye!")
            sys.exit(0)

        if choice in ("", "1"):
            while True:
                pi_res, _ = run_single_calculation(
                    digits=digits,
                    algorithm=algorithm,
                    use_code_map=use_code_map,
                    delay=delay,
                )
                last_pi_str = pi_res
                next_action = post_calc_actions(
                    digits=digits,
                    algorithm=algorithm,
                    use_code_map=use_code_map,
                    delay=delay,
                    pi_str=last_pi_str,
                )
                if next_action == "rerun":
                    continue
                elif next_action == "menu":
                    break
                elif next_action == "quit":
                    clear_screen()
                    print("Session terminated. Goodbye!")
                    sys.exit(0)

        elif choice == "2":
            clear_screen()
            print("=" * 66)
            print("                     TARGET PRECISION                         ")
            print("=" * 66)
            print("  [1] 1,000 digits    (Instant: < 0.01s)")
            print("  [2] 10,000 digits   (Fast: ~0.02s)")
            print("  [3] 50,000 digits   (High: ~0.11s)")
            print("  [4] 100,000 digits  (Deep: ~0.31s)")
            print("  [5] 250,000 digits  (Ultra: ~1.4s)")
            print("  [C] Custom digit count")
            print("  [B] Back to Main Menu")
            print("-" * 66)
            p_choice = input("Choice [1-5, C, or B]: ").strip().lower()
            presets = {"1": 1000, "2": 10000, "3": 50000, "4": 100000, "5": 250000}
            if p_choice in presets:
                digits = presets[p_choice]
                status_msg = f"\033[92m✔ Target precision set to {digits:,} digits.\033[0m"
            elif p_choice == "c":
                cust = input("Enter custom number of digits: ").strip()
                try:
                    val = int(cust)
                    if val > 0:
                        digits = val
                        status_msg = f"\033[92m✔ Target precision set to {digits:,} digits.\033[0m"
                    else:
                        status_msg = "\033[93m⚠ Invalid number.\033[0m"
                except ValueError:
                    status_msg = "\033[93m⚠ Invalid input.\033[0m"

        elif choice == "3":
            clear_screen()
            print("=" * 66)
            print("                    DECODING ALGORITHM                        ")
            print("=" * 66)
            print("  [1] Chudnovsky (BS)  (Integer Binary Splitting) [Fastest]")
            print("  [2] Gauss-Legendre   (Salamin-Brent AGM, dynamic precision)")
            print("  [3] Gregory-Leibniz  (Paired alternating series)")
            print("  [B] Back to Main Menu")
            print("-" * 66)
            a_choice = input("Choice [1-3, or B]: ").strip()
            if a_choice == "1":
                algorithm = "chudnovsky"
                status_msg = f"\033[92m✔ Switched algorithm to Chudnovsky (Binary Splitting).\033[0m"
            elif a_choice == "2":
                algorithm = "gauss-legendre"
                status_msg = f"\033[92m✔ Switched algorithm to Gauss-Legendre.\033[0m"
            elif a_choice == "3":
                algorithm = "leibniz"
                status_msg = f"\033[92m✔ Switched algorithm to Gregory-Leibniz.\033[0m"

        elif choice == "4":
            clear_screen()
            print("=" * 66)
            print("               CODE MAP & PACING CONFIGURATION                ")
            print("=" * 66)
            print(f"  [1] Toggle Code Map on/off       (Currently: {map_status})")
            print(f"  [2] Toggle Live Pacing on/off    (Currently: {'ON' if delay > 0 else 'OFF'})")
            print("  [3] Fast Pacing Preset      (0.010s / frame)")
            print("  [4] Normal Pacing Preset    (0.030s / frame)")
            print("  [5] Slow Pacing Preset      (0.080s / frame)")
            print("  [6] Enter Custom Delay in seconds")
            print("  [B] Back to Main Menu")
            print("-" * 66)
            s_choice = input("Choice [1-6, or B]: ").strip()
            if s_choice == "1":
                use_code_map = not use_code_map
                status_msg = f"\033[92m✔ Code Map {'ENABLED' if use_code_map else 'DISABLED'}.\033[0m"
            elif s_choice == "2":
                delay = 0.02 if delay == 0.0 else 0.0
                status_msg = f"\033[92m✔ Pacing {'ENABLED (0.020s)' if delay > 0 else 'DISABLED'}.\033[0m"
            elif s_choice == "3":
                delay = 0.01
                status_msg = "\033[92m✔ Pacing set to 0.010s.\033[0m"
            elif s_choice == "4":
                delay = 0.03
                status_msg = "\033[92m✔ Pacing set to 0.030s.\033[0m"
            elif s_choice == "5":
                delay = 0.08
                status_msg = "\033[92m✔ Pacing set to 0.080s.\033[0m"
            elif s_choice == "6":
                c_delay = input("Enter delay per frame in seconds (e.g. 0.02): ").strip()
                try:
                    delay = max(0.0, float(c_delay))
                    status_msg = f"\033[92m✔ Pacing set to {delay:.3f}s.\033[0m"
                except ValueError:
                    status_msg = "\033[93m⚠ Invalid input.\033[0m"

        elif choice == "5":
            benchmark_algorithms(digits=min(digits, 10000))
            input("\nPress Enter to return to menu...")

        elif choice == "6":
            if not last_pi_str:
                clear_screen()
                print("Decoding Pi stream first...")
                last_pi_str, _ = run_single_calculation(
                    digits=digits,
                    algorithm=algorithm,
                    use_code_map=use_code_map,
                    delay=0.0,
                )
            search_in_pi(last_pi_str)
            input("\nPress Enter to return to menu...")

        elif choice == "7":
            clear_screen()
            print("=" * 66)
            print("                  SAVE STREAM TO FILE                         ")
            print("=" * 66)
            if not last_pi_str:
                print("No calculated data in buffer. Run a decode first.")
            else:
                fname = input(f"Filename [default: pi_{digits}_digits.txt]: ").strip()
                if not fname:
                    fname = f"pi_{digits}_digits.txt"
                with open(fname, "w", encoding="utf-8") as f:
                    f.write(last_pi_str)
                status_msg = f"\033[92m✔ Stream exported to {fname}\033[0m"
            if not status_msg:
                input("\nPress Enter to return to menu...")

        elif choice in ("0", "q", "exit"):
            clear_screen()
            print("Session terminated. Goodbye!")
            sys.exit(0)
        else:
            status_msg = "\033[93m⚠ Invalid choice. Select 0-7.\033[0m"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="High-performance Pi (π) calculator with Orbital Code Map telemetry. Created by Dreamlexxi."
    )
    parser.add_argument(
        "digits",
        nargs="?",
        type=int,
        default=None,
        help="Number of decimal digits to compute (e.g. 10000, 100000).",
    )
    parser.add_argument(
        "-a",
        "--algorithm",
        choices=["chudnovsky", "gauss-legendre", "leibniz"],
        default="chudnovsky",
        help="Algorithm (default: chudnovsky).",
    )
    parser.add_argument(
        "--menu",
        action="store_true",
        help="Open interactive console dashboard.",
    )
    parser.add_argument(
        "--animate",
        action="store_true",
        help="Enable pacing animation for the code map display.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.0,
        help="Custom delay per frame in seconds (e.g. 0.02).",
    )
    parser.add_argument(
        "--no-map",
        action="store_true",
        help="Disable the orbital code map progress display.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Save calculated digits to file.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.menu or (args.digits is None and sys.stdin.isatty()):
        interactive_menu()
        return

    digits = args.digits if args.digits is not None else 10000
    if digits < 1:
        print("Error: Digits must be an integer greater than 0.")
        sys.exit(1)

    delay = args.delay
    if args.animate and delay == 0.0:
        delay = 0.03

    pi_str, _ = run_single_calculation(
        digits=digits,
        algorithm=args.algorithm,
        use_code_map=not args.no_map,
        delay=delay,
        output_file=args.output,
    )

    if sys.stdin.isatty():
        action = post_calc_actions(
            digits=digits,
            algorithm=args.algorithm,
            use_code_map=not args.no_map,
            delay=delay,
            pi_str=pi_str,
        )
        if action == "rerun":
            while True:
                p_str, _ = run_single_calculation(
                    digits=digits,
                    algorithm=args.algorithm,
                    use_code_map=not args.no_map,
                    delay=delay,
                    output_file=args.output,
                )
                act = post_calc_actions(
                    digits=digits,
                    algorithm=args.algorithm,
                    use_code_map=not args.no_map,
                    delay=delay,
                    pi_str=p_str,
                )
                if act == "rerun":
                    continue
                elif act == "menu":
                    interactive_menu()
                    break
                elif act == "quit":
                    clear_screen()
                    break
        elif action == "menu":
            interactive_menu()


if __name__ == "__main__":
    main()
