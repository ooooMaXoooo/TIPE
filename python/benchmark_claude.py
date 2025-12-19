import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import TIPE_SimuOrbit
import os
import json
import csv
from datetime import datetime
from dataclasses import dataclass, asdict
from typing import List, Tuple, Dict, Any

# Forcer le backend interactif QtAgg
matplotlib.use("QtAgg")

@dataclass
class BenchmarkResult:
    """Classe pour stocker les résultats d'un benchmark"""
    n_points: int
    n_orbits: int
    time_cpp: float
    time_conversion: float
    time_numpy: float
    time_total_cpp: float
    throughput_cpp: float
    throughput_numpy: float
    speedup: float
    memory_usage_mb: float = 0.0

class BenchmarkRunner:
    """Classe principale pour gérer les benchmarks"""
    
    def __init__(self, mu: float = 1.32712440018e20, R_Mars: float = 2.28e11):
        self.mu = mu
        self.R_Mars = R_Mars
        self.r1 = [1.5e11, 0.0, 0.0]
        self.results: List[BenchmarkResult] = []
        self.benchmark_dir = None
        self.timestamp = None
        
    def setup_benchmark_folder(self) -> Tuple[str, str]:
        """Création du dossier benchmark avec timestamp"""
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.benchmark_dir = f"benchmarks/benchmark_{self.timestamp}"
        os.makedirs(self.benchmark_dir, exist_ok=True)
        print(f"📁 Dossier de sauvegarde: {self.benchmark_dir}/")
        print(f"🕐 Timestamp: {self.timestamp}")
        return self.benchmark_dir, self.timestamp

    def prepare_data(self, n_points: int) -> Tuple[List[List[float]], List[float]]:
        """Préparation des données d'entrée pour le benchmark"""
        tof_list = np.linspace(20, 500, n_points)
        theta_list = np.linspace(0, 2*np.pi, n_points)
        Theta, ToF = np.meshgrid(theta_list, tof_list)
        
        r2_list = []
        tof_flat = []
        
        for i in range(n_points):
            for j in range(n_points):
                theta_rad = Theta[i, j]
                r2 = self.R_Mars * np.array([np.cos(theta_rad), np.sin(theta_rad), 0.0])
                r2_list.append(r2.tolist())
                tof_flat.append(tof_list[i])
                
        return r2_list, tof_flat

    def estimate_memory_usage(self, n_points: int) -> float:
        """Estimation de l'usage mémoire en MB"""
        n_orbits = n_points * n_points
        # Estimation approximative : 3 vecteurs 3D + 1 scalaire par orbite
        bytes_per_orbit = 3 * 3 * 8 + 8  # 8 bytes par float64
        total_bytes = n_orbits * bytes_per_orbit * 2  # facteur 2 pour les structures temporaires
        return total_bytes / (1024 * 1024)

    def run_single_benchmark(self, n_points: int) -> BenchmarkResult:
        """Exécute un benchmark pour une taille de grille donnée"""
        print(f"\n=== Benchmark pour N={n_points} ===")
        r2_list, tof_flat = self.prepare_data(n_points)
        n_orbits = n_points * n_points

        # Estimation mémoire
        memory_estimate = self.estimate_memory_usage(n_points)
        print(f"  📊 Orbites à calculer: {n_orbits:,}")
        print(f"  💾 Mémoire estimée: {memory_estimate:.1f} MB")

        # lambert_batch (C++ → list)
        print("  🔄 Calcul lambert_batch (C++ → list)...")
        start_time = time.perf_counter()
        res_vec_list = TIPE_SimuOrbit.orbit.lambert_batch(self.r1, r2_list, tof_flat, self.mu)
        time_cpp = time.perf_counter() - start_time

        # Conversion en numpy
        print("  🔄 Conversion vers tableau numpy...")
        start_time = time.perf_counter()
        res_vec = np.array(res_vec_list)
        time_conversion = time.perf_counter() - start_time

        time_total_cpp = time_cpp + time_conversion

        # lambert_batch_numpy
        print("  🔄 Calcul lambert_batch_numpy...")
        start_time = time.perf_counter()
        _ = TIPE_SimuOrbit.orbit.lambert_batch_numpy(self.r1, r2_list, tof_flat, self.mu)
        time_numpy = time.perf_counter() - start_time

        # Calculs des métriques
        throughput_cpp = n_orbits / time_total_cpp
        throughput_numpy = n_orbits / time_numpy
        speedup = time_total_cpp / time_numpy
        numpy_advantage = time_numpy / time_total_cpp  # < 1 = NumPy plus rapide

        result = BenchmarkResult(
            n_points=n_points,
            n_orbits=n_orbits,
            time_cpp=time_cpp,
            time_conversion=time_conversion,
            time_numpy=time_numpy,
            time_total_cpp=time_total_cpp,
            throughput_cpp=throughput_cpp,
            throughput_numpy=throughput_numpy,
            speedup=speedup,
            memory_usage_mb=memory_estimate
        )

        print(f"  ✅ C++ + conversion: {time_total_cpp:.3f}s ({throughput_cpp:,.0f} orbites/s)")
        print(f"  ✅ NumPy direct: {time_numpy:.3f}s ({throughput_numpy:,.0f} orbites/s)")
        print(f"  📈 Facteur NumPy: x{numpy_advantage:.3f} ({'NumPy plus rapide' if numpy_advantage < 1 else 'C++ plus rapide'})")

        return result

    def run_full_benchmark(self, grid_sizes: List[int]) -> List[BenchmarkResult]:
        """Exécute le benchmark complet pour toutes les tailles"""
        print("🚀 DÉBUT DU BENCHMARK")
        print("="*60)
        
        self.setup_benchmark_folder()
        self.results = []

        for n in grid_sizes:
            try:
                result = self.run_single_benchmark(n)
                self.results.append(result)
            except Exception as e:
                print(f"❌ Erreur pour N={n}: {e}")
                continue

        print(f"\n🎯 BENCHMARK TERMINÉ - {len(self.results)} tests réussis")
        print("="*60)
        return self.results

    def save_results_json(self, filename: str = None):
        """Sauvegarde les résultats en JSON"""
        if filename is None:
            filename = f"results_{self.timestamp}.json"
        filepath = os.path.join(self.benchmark_dir, filename)
        
        data = {
            'metadata': {
                'timestamp': self.timestamp,
                'mu': self.mu,
                'R_Mars': self.R_Mars,
                'r1': self.r1,
                'total_tests': len(self.results)
            },
            'results': [asdict(result) for result in self.results]
        }
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        print(f"💾 Résultats sauvegardés: {filepath}")

    def save_results_csv(self, filename: str = None):
        """Sauvegarde les résultats en CSV"""
        if filename is None:
            filename = f"results_{self.timestamp}.csv"
        filepath = os.path.join(self.benchmark_dir, filename)
        
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            if self.results:
                writer = csv.DictWriter(f, fieldnames=asdict(self.results[0]).keys())
                writer.writeheader()
                for result in self.results:
                    writer.writerow(asdict(result))
        print(f"📊 CSV sauvegardé: {filepath}")

    def print_summary_table(self):
        """Affiche un tableau de synthèse des résultats"""
        if not self.results:
            print("❌ Aucun résultat à afficher")
            return

        print("\n📋 RÉSUMÉ DES PERFORMANCES")
        print("="*90)
        print(f"{'N':>5} {'Orbites':>10} {'C++ (s)':>10} {'Conv (s)':>10} {'NumPy (s)':>10} {'Facteur NumPy':>15} {'Débit NumPy':>15}")
        print("-"*90)
        
        for result in self.results:
            numpy_factor = result.time_numpy / result.time_total_cpp
            factor_str = f"x{numpy_factor:.3f}"
            
            print(f"{result.n_points:>5} "
                  f"{result.n_orbits:>10,} "
                  f"{result.time_cpp:>10.3f} "
                  f"{result.time_conversion:>10.3f} "
                  f"{result.time_numpy:>10.3f} "
                  f"{factor_str:>15} "
                  f"{result.throughput_numpy:>15,.0f}")

    def create_enhanced_plots(self):
        """Création des graphiques améliorés"""
        if not self.results:
            print("❌ Aucun résultat pour les graphiques")
            return None

        # Extraction des données
        grid_sizes = [r.n_points for r in self.results]
        times_cpp = [r.time_cpp for r in self.results]
        times_conv = [r.time_conversion for r in self.results]
        times_numpy = [r.time_numpy for r in self.results]
        times_total_cpp = [r.time_total_cpp for r in self.results]
        throughput_cpp = [r.throughput_cpp for r in self.results]
        throughput_numpy = [r.throughput_numpy for r in self.results]

        # Configuration du style et des positions
        try:
            plt.style.use('seaborn-v0_8-darkgrid')
        except:
            plt.style.use('default')
            
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        plt.subplots_adjust(hspace=0.3, wspace=0.3)
        
        # Positions communes pour tous les graphiques avec barres
        x_positions = np.arange(len(grid_sizes))
        bar_width = 0.6
        
        # Couleurs cohérentes
        colors = {
            'cpp': '#FF6B6B',
            'conv': '#4ECDC4', 
            'numpy': '#45B7D1',
            'total': '#96CEB4'
        }

        # 1) Comparaison directe des temps (barres côte à côte au lieu d'empilées)
        bar_width_comp = 0.35
        x_cpp = x_positions - bar_width_comp/2
        x_numpy = x_positions + bar_width_comp/2
        
        bars_cpp = axes[0,0].bar(x_cpp, times_total_cpp, bar_width_comp, 
                                 label="C++ + conversion", color=colors['cpp'], alpha=0.8)
        bars_numpy = axes[0,0].bar(x_numpy, times_numpy, bar_width_comp, 
                                   label="NumPy direct", color=colors['numpy'], alpha=0.8)
        
        # Annotations avec gain en % sur les barres NumPy (toujours affichées)
        for i, (cpp_time, numpy_time) in enumerate(zip(times_total_cpp, times_numpy)):
            gain_percent = (cpp_time - numpy_time) / cpp_time * 100
            # Afficher toutes les annotations (pas de seuil)
            color = 'green' if gain_percent > 0 else 'red'
            axes[0,0].annotate(f'{gain_percent:+.1f}%', 
                              xy=(x_numpy[i], numpy_time), 
                              xytext=(0, 8), textcoords='offset points',
                              ha='center', va='bottom', fontweight='bold',
                              color=color, fontsize=9)
        
        axes[0,0].set_xticks(x_positions)
        axes[0,0].set_xticklabels([str(n) for n in grid_sizes])
        axes[0,0].set_xlabel("Taille de grille N", fontsize=12)
        axes[0,0].set_ylabel("Temps (s)", fontsize=12)
        axes[0,0].set_title("Comparaison directe des temps d'exécution", fontsize=14, fontweight='bold')
        axes[0,0].legend()
        axes[0,0].grid(axis='y', alpha=0.3)

        # 2) Facteur NumPy vs C++ AVEC MOYENNE
        numpy_factors = [numpy/total for numpy, total in zip(times_numpy, times_total_cpp)]
        
        # Ligne de référence à 1.0 (performance égale)
        axes[0,1].axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Performance égale')
        
        # Ligne de moyenne des performances
        mean_factor = np.mean(numpy_factors)
        axes[0,1].axhline(y=mean_factor, color='red', linestyle=':', alpha=0.8, linewidth=2,
                         label=f'Moyenne: x{mean_factor:.3f}')
        
        # Courbe des données
        axes[0,1].plot(grid_sizes, numpy_factors, "o-", linewidth=3, markersize=8, 
                      color=colors['numpy'], label='Facteur NumPy')
        
        # Annotations des valeurs - FORMAT x0.XXX
        for x, y in zip(grid_sizes, numpy_factors):
            axes[0,1].annotate(f'x{y:.3f}', xy=(x, y), xytext=(5, 10), 
                              textcoords='offset points', ha='left', fontweight='bold')
        
        axes[0,1].set_xlabel("Taille de grille N", fontsize=12)
        axes[0,1].set_ylabel("Facteur NumPy/C++", fontsize=12)
        axes[0,1].set_title("Performance relative NumPy vs C++", fontsize=14, fontweight='bold')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)

        # 3) Débit avec échelle log
        axes[1,0].loglog(grid_sizes, throughput_cpp, "s-", linewidth=2, markersize=8, 
                        color=colors['cpp'], label="C++ + conversion")
        axes[1,0].loglog(grid_sizes, throughput_numpy, "o-", linewidth=2, markersize=8, 
                        color=colors['numpy'], label="NumPy direct")
        
        axes[1,0].set_xlabel("Taille de grille N", fontsize=12)
        axes[1,0].set_ylabel("Débit (orbites/s)", fontsize=12)
        axes[1,0].set_title("Débit de calcul (échelle log)", fontsize=14, fontweight='bold')
        axes[1,0].legend()
        axes[1,0].grid(True, alpha=0.3)

        # 4) Barres avec facteurs NumPy ET DEBUG
        bars = axes[1,1].bar(x_positions, numpy_factors, width=bar_width, 
                            color=colors['numpy'], alpha=0.8)
        
        # Ligne de référence à 1.0 (performance égale)
        axes[1,1].axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='Performance égale')
        
        # DEBUG: afficher les valeurs des barres
        print("\n🔍 DEBUG - Valeurs des barres (graphique 4):")
        for i, (factor, bar) in enumerate(zip(numpy_factors, bars)):
            height = bar.get_height()
            print(f"  Barre {i}: factor={factor:.4f}, height={height:.4f}")
        
        # Annotations DANS les barres - FORMAT SIMPLE x0.XXX
        for i, (bar, factor) in enumerate(zip(bars, numpy_factors)):
            # Récupération des dimensions de la barre
            height = bar.get_height()
            x_center = bar.get_x() + bar.get_width() / 2
            y_center = height / 2  # Centre vertical de la barre
            
            print(f"  Position texte {i}: x={x_center:.3f}, y={y_center:.3f}")
            
            # Texte à afficher - FORMAT SIMPLE x0.XXX
            text = f'x{factor:.3f}'
            
            # Placement du texte au centre de la barre
            axes[1,1].text(x_center, y_center, text, 
                          ha='center', va='center', 
                          fontsize=10, fontweight='bold',
                          color='black',
                          bbox=dict(boxstyle="round,pad=0.3", 
                                   facecolor='white', alpha=0.9, 
                                   edgecolor='black', linewidth=1))
        
        axes[1,1].set_xticks(x_positions)
        axes[1,1].set_xticklabels([str(n) for n in grid_sizes])
        axes[1,1].set_xlabel("Taille de grille N", fontsize=12)
        axes[1,1].set_ylabel("Facteur NumPy/C++", fontsize=12)
        axes[1,1].set_title("Facteur de performance NumPy", fontsize=14, fontweight='bold')
        axes[1,1].legend()
        axes[1,1].grid(axis='y', alpha=0.3)

        # Sauvegarde avec titre général
        plt.suptitle(f'Benchmark Orbites - {self.timestamp}', fontsize=16, fontweight='bold', y=0.98)
        
        png_file = os.path.join(self.benchmark_dir, f"benchmark_{self.timestamp}.png")
        pdf_file = os.path.join(self.benchmark_dir, f"benchmark_{self.timestamp}.pdf")
        
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.savefig(pdf_file, bbox_inches='tight')
        print(f"📊 Graphiques sauvegardés: {png_file}")
        print(f"📊 PDF sauvegardé: {pdf_file}")

        return fig

def main():
    """Fonction principale"""
    # Configuration du benchmark
    benchmark = BenchmarkRunner()
    
    # Tailles de grille à tester
    grid_sizes = [50*i for i in range(1, 11)]
    
    # Exécution du benchmark
    results = benchmark.run_full_benchmark(grid_sizes)
    
    if results:
        # Affichage du résumé
        benchmark.print_summary_table()
        
        # Sauvegarde des données
        benchmark.save_results_json()
        benchmark.save_results_csv()
        
        # Création et affichage des graphiques
        try:
            fig = benchmark.create_enhanced_plots()
            if fig:
                plt.show()
        except Exception as e:
            print(f"❌ Erreur lors de la création des graphiques: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("❌ Aucun résultat à traiter")

if __name__ == "__main__":
    main()