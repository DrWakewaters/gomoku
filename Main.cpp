#include "Board.h"
#include "Game.h"
#include "Player.h"
#include "AI.h"

// Compile with -std=c++11 -Ofast.
int main() {
	int8_t maxSearchDepth = 16;
	int8_t boardHeight = 19;
	int8_t boardWidth = 19;
	// This is not actually the max time for an AI player to make a ply. Rather, it can start another iteration of the iterative deepener when the computation time is smaller than this, so the total time could, in rare, extreme cases, be 100 timer longer (typically it is less than 10 times longer, and often just 1-5 times longer).
	double maxSearchTime = 1.0;
	bool blackIsHuman = true;
	bool whiteIsHuman = false;
	bool printAIInformation = true;
	Game game(maxSearchDepth, boardHeight, boardWidth, maxSearchTime, blackIsHuman, whiteIsHuman, printAIInformation);
	bool blackWon = false;
	bool whiteWon = false;
	std::cout << "When making a ply, write it as \"y x\", without the quotation marks." << std::endl;
	game.printBoard();
	while(true) {
		if(game.boardIsFull()) {
			std::cout << "It is a draw." << std::endl;
			break;
		}
		blackWon = game.makeBlackPly();
		AI::evaluateBoard(game.board, false, true);
		if(blackWon) {
			std::cout << "Black won after " << game.getNumberOfPliesMade() << " plies." << std::endl;
			break;
		}
		if(game.boardIsFull()) {
			std::cout << "It is a draw." << std::endl;
			break;
		}
		whiteWon = game.makeWhitePly();
		AI::evaluateBoard(game.board, true, true);
		if(whiteWon) {
			std::cout << "White won after " << game.getNumberOfPliesMade() << " plies." << std::endl;
			break;
		}
	}
	return 0;
}
