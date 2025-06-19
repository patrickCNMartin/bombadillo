{
  description = "Bombadillo";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }: let
    supportedSystems = [ "aarch64-darwin" "x86_64-linux" ];
    forAllSystems = f: nixpkgs.lib.genAttrs supportedSystems (system: f { inherit system; });
  in {
    devShells = forAllSystems ({ system }: let
      pkgs = import nixpkgs { inherit system; };
    in {
      default = pkgs.mkShell {
        buildInputs = [
          pkgs.julia-bin
        ];
      };
    });
  };
}
